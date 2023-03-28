#include "Model_LocalBranching.h"

#define EPS 1.0e-6

ILOSTLBEGIN
//typedef IloArray<IloNumArray>    Num2D;
//typedef IloArray<IloNumVarArray>    NumVar2D;
//typedef IloArray<NumVar2D>    NumVar3D;
//typedef IloArray<IloBoolVarArray> BoolVar2D;
//typedef IloArray<BoolVar2D> BoolVar3D;


void Model_LocalBranching::SetData(const Data_Improved& newData)
{
	data = newData;
}


void Model_LocalBranching::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);
	TimeLimit = params.value("timelimit", 7200);


	//Define constants for easier reading and writing
	T = data.T;
	M_bar = data.M_bar;
	P = data.P;

	//Initialise data-dependent structures
	BestSolution = ArrayXXd::Constant(T, M_bar, 0.0);
	core_point = ArrayXXd::Constant(T, M_bar, 0.0);
	for (int t = 0; t < T; t++) {
		core_coverage.push_back(VectorXd::Zero(P[t]));
	}
	


	env = IloEnv();
	model = IloModel(env);
	IloCplex cplex(model);


	//Set parameters
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::TimeLimit, TimeLimit);

	cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
	cplex.setParam(IloCplex::Param::WorkMem, 100000);


	cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
	cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 1);

	cplex.setParam(IloCplex::Param::Emphasis::MIP, 2);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	//cplex.setParam(IloCplex::Param::Emphasis::MIP, 0);

	cplex.setParam(IloCplex::Param::ParamDisplay, 0);
	cplex.setParam(IloCplex::Param::Read::WarningLimit, 0);
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	//cplex.setParam(IloCplex::Param::MIP::Interval, 10000);
	cplex.setParam(IloCplex::Param::Simplex::Display, 2);
		

	//Create variables
	////Charging station placement (where j_bar defines a pair (j,k) in data)
	x = BoolVar2D(env, T);
	for (int t = 0; t < T; t++) {
		x[t] = IloBoolVarArray(env, M_bar);
	}

	////Auxiliary variable for Bender's cuts, given for each time period
	theta= IloNumVarArray(env, T, 0.0, IloInfinity);

	////Auxiliary variable for overall objective
	theta_obj = IloNumVar(env, 0.0, IloInfinity);

	////Pregenerated array of binary variables for use in generating trust cuts
	////This uses a fixed size, to avoid stale reference potential when modifying the size
	trust = BoolVar2D(env, 1024);
	for (int i = 0; i < 1024; i++) {
		trust[i] = IloBoolVarArray(env, T);
	}


	//Constraints
	////Budget, year 0
	IloExpr budget0(env);
	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		int j = data.params["station_coord"][j_bar][0];
		int k = data.params["station_coord"][j_bar][1];
		budget0 += (double)data.params["c"][0][j][k] * (x[0][j_bar] - (int)data.params["x0"][j][k]);
	}
	model.add(budget0 <= (double)data.params["B"][0]);
	budget0.end();

	////Budget, year 1+
	for (int t = 1; t < T; t++) {
		IloExpr budget(env);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			budget += (double)data.params["c"][t][j][k] * (x[t][j_bar] - x[t - 1][j_bar]);
		}
		model.add(budget <= (double)data.params["B"][t]);
		budget.end();
	}

	////Can't remove outlets, year 0
	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		int j = data.params["station_coord"][j_bar][0];
		int k = data.params["station_coord"][j_bar][1];
		model.add(x[0][j_bar] >= (int)data.params["x0"][j][k]);
	}

	////Can't remove outlets, year 1+
	for (int t = 1; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			model.add(x[t][j_bar] >= x[t - 1][j_bar]);
		}
	}

	////At least k outlets, year 0+
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int k = data.params["station_coord"][j_bar][1];
			if (k > 1) {
				model.add(x[t][j_bar] <= x[t][j_bar - 1]);
			}
		}
	}

	///Redundant, but prevents crash from unextracted variables
	model.add(theta_obj >= 0);
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			model.add(x[t][j_bar] >= 0);
		}
		model.add(theta[t] >= 0);

		for (int i = 0; i < 1024; i++) {
			model.add(trust[i][t] >= 0);
		}
	}


	//Create warmstart solution (greedy)
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	vector<vector<int>> sol;
	Greedy_Improved G;
	G.SetData(data);
	G.Solve(false);
	sol = G.Solution;
	LowerBound = G.SolutionQuality;

	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			int val = 0;
			if (sol[t][j] >= k) {
				val = 1;
			}
			startVar.add(x[t][j_bar]);
			startVal.add(val);
			BestSolution(t, j_bar) = val;
		}

		startVar.add(theta[t]);
		IloNum val = 0;
		VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			if (sol[t][j] >= k) {
				I_tilde += data.a[t].col(j_bar);
			}
		}
		val += data.weights[t].dot((VectorXd)(I_tilde.array() >= 1).matrix().cast<double>());
		startVal.add(val);
	}


	cplex.addMIPStart(startVar, startVal);
	startVal.end();
	startVar.end();



	////Upper bound for theta (ensures bounded problem)
	for (int t = 0; t < T; t++) {
		model.add(theta[t] <= data.weights[t].sum());
	}


	//Objective
	if (verbose) { cplex.out() << "Adding objective \n"; }
	IloExpr obj(env);
	for (int t = 0; t < T; t++) {
		obj += theta[t] + data.Precovered[t];
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			obj += data.Ps[t](j_bar) * x[t][j_bar];
		}
	}
	model.add(theta_obj == obj);
	obj.end();

	//Set objective via proxy
	model.add(IloMaximize(env, theta_obj));

	//Calculate solution of LP relaxation and add associated fractional Bender's cuts to model
	{
		vector<IloConversion> conv_x;
		for (int t = 0; t < T; t++) {
			IloConversion conv(env, x[t], IloNumVar::Type::Float);
			model.add(conv);
			conv_x.push_back(conv);
		}
			
		Callback_Improved cb(data, x, theta);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);
		time_t start;
		time(&start);
		cplex.solve();
		TimeLimit -= time(NULL) - start;

		UpperBound = cplex.getBestObjValue();
		
		for (auto conv : conv_x) {
			model.remove(conv);
		}
	}

	//Link callback
	Callback_LocalBranching cb(data, cplex, x, theta); 
	cb.thresholdObjective = LowerBound + 1e-4 * (LowerBound + 1e-10);
	cb.LowerBound = LowerBound;
		
	CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
		| IloCplex::Callback::Context::Id::Relaxation
		| IloCplex::Callback::Context::Id::Branching;
	cplex.use(&cb, contextmask);


	
	//Solve and get results
	time_t start;
	time(&start);
	while (SolveTime < TimeLimit) {
		try {
			bool solved = cplex.solve();
		}
		catch (IloException & e) {
			cplex.out() << "Exception: " << e << endl;
		}
		if (cb.status == CallbackStatus::Optimal) {
			cout << "CPLEX optimality conditions reached." << endl;
			break;
		}
		UpperBound = cb.UpperBound;
		IntegerCutCallback(cb.x_tilde, cb.theta_tilde, -1, true, cb.status);
		if (UpdatedBound) {
			cb.thresholdObjective = thresholdObjective;
			cb.LowerBound = LowerBound;
			UpdatedBound = false;
		}

		cout << "Upper bound: " << UpperBound << "   Threshold: " << thresholdObjective << endl;
		if (UpperBound < thresholdObjective) { cout << "Optimality gap reached." << endl;  break; }

		cb.status = CallbackStatus::Unknown;
		SolveTime = time(NULL) - start;
		cout << "Current time: " << SolveTime << endl;
		if (7200 - SolveTime > 1) { cplex.setParam(IloCplex::Param::TimeLimit, TimeLimit - SolveTime); }
	}
		


	if (verbose) {
		cplex.out() << "Solve time: " << SolveTime << endl;
		cplex.out() << "Objective value: " << LowerBound << endl;
		cplex.out() << "Upper bound: " << UpperBound << endl;
		cplex.out() << "Optimality gap: " << (UpperBound - LowerBound)/LowerBound << endl;
		cplex.out() << "\n" << endl;
	}


	IloAlgorithm::Status status = cplex.getStatus();
	stats["CplexStatus"] = status;
	stats["SolveTime (x100)"] = (int)100 * SolveTime;
	stats["ObjectiveValue (x100)"] = (int)100 * LowerBound;


	cplex.clearUserCuts();
	cplex.end();
	model.end();
	env.end();
	//fout.close();

}



void Model_LocalBranching::GetSolution(IloCplex& cplex, BoolVar2D& x)
{
	for (int t = 0; t < data.T; t++) {
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			BestSolution(t, j_bar) = cplex.getValue(x[t][j_bar]);
		}
	}
}



//Updates the core point used in Bender's cuts and the associated coverage
void Model_LocalBranching::UpdateCorePoint(const ArrayXXd& x_new)
{
	core_point = (core_point + x_new) / 2;
	for (int t = 0; t < T; t++) {
		VectorXd cover = VectorXd::Constant(P[t], 0.0);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (core_point(t, j_bar) > EPS) {
				cover += core_point(t, j_bar) * data.a[t].col(j_bar);
			}
		}
		core_coverage[t] = cover;
	}
}


//Add Bender's cuts to model (the one saved as class attribute), either using Multi1B1 or Multi1PO1 cuts
void Model_LocalBranching::AddBendersCuts(const ArrayXXd& x_tilde, const ArrayXd& theta_tilde)
{
	vector<VectorXd> I_tilde = CalculateItilde(x_tilde);

	////Multi-cut (by year), B1
//for (int t = 0; t < T; t++) {
//	IloExpr lhs(env);
//	lhs -= theta[t];
//	double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
//	VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
//	}

//	switch (context.getId()) {
//	case (IloCplex::Callback::Context::Id::Candidate):
//		if (context.getCandidateValue(lhs) < -covered -EPS) {
//			//nCutsAdded += 1;
//			context.rejectCandidate(lhs >= -covered);
//		}
//		lhs.end();
//		break;
//	case (IloCplex::Callback::Context::Id::Relaxation):
//		if (context.getRelaxationValue(lhs) < -covered -EPS) {
//			//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
//			context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
//		}
//		lhs.end();
//		break;
//	}
//}



//Multi-cut (by year), Multi1PO1
	for (int t = 0; t < T; t++) {
		double eval = 0.0;
		IloExpr lhs(env);
		lhs -= theta[t];
		eval -= theta_tilde(t);
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

		eval += covered;

		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			double coef = (data.CutCoeffs[t].col(j_bar).dot(uncovered));
			lhs +=  coef * x[t][j_bar];
			if (x_tilde(t, j_bar) > 1 - EPS) { eval += coef; }
		}
		if (eval < -EPS) {
			model.add(lhs >= -covered);
		}
	}

}


//Add trust cuts to *submodel* (not necessarily the class attribute), imposing that we must be within *distance* from *x_tilde*
void Model_LocalBranching::AddTrustCutsInner(IloModel& submodel, const ArrayXXd& x_tilde, double distance)
{
	for (int t = 0; t < T; t++) {
		IloExpr delta(env);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
			else { delta += x[t][j_bar]; }
		}
		submodel.add(delta <= distance);
		delta.end();
	}
}


//Add trust cuts to model (the one saved as class attribute), imposing that we must be at least *distance* aways from *x_tilde*
void Model_LocalBranching::AddTrustCutsOuter(const ArrayXXd& x_tilde, double distance)
{
	//Trust count should be gated by a lock, but in single-threaded mode this is not necessary
	trustCount += 1;
	if (trustCount >= trust.getSize() - 1) {
		cout << "Maximum number of trust constraints exceeded." << endl;
	}

	for (int t = 0; t < T; t++) {
		IloExpr delta(env);
		delta += distance * trust[trustCount][t];
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
			else { delta += x[t][j_bar]; }
		}
		model.add(delta >= distance);
		delta.end();
	}

	model.add(IloSum(trust[trustCount]) <= T - 1); //This constraint ensures that at least one of the distances is >= *distance*
}


//Add simple taboo (/no good) cut, imposing that we may not select *x_tilde* again
void Model_LocalBranching::AddTabooCut(IloModel& submodel, const ArrayXXd& x_tilde)
{
	IloExpr taboo(env);
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { taboo += 1 - x[t][j_bar]; }
			else { taboo += x[t][j_bar]; }
		}
	}
	submodel.add(taboo >= 1);
}



void Model_LocalBranching::IntegerCutCallback(ArrayXXd& x_tilde, ArrayXd& theta_tilde, double obj_tilde, bool diversify, CallbackStatus status)
{

	ArrayXXd x_hat = ArrayXXd::Zero(T, M_bar);
	if (obj_tilde == -1) { obj_tilde = CalculateObjective(x_tilde); }
	if ((status != CallbackStatus::FractionalRound) && (obj_tilde > LowerBound)) { cout << "Lower bound improved: " << obj_tilde << endl;  UpdateSolution(x_tilde, obj_tilde); }
	double obj_hat = -1.0;

	bool foundImprovement = false;
	do
	{
		foundImprovement = false;
		IloAlgorithm::Status status = SubproblemRestricted(x_tilde, x_hat, obj_hat);
		switch (status)
		{
			////////////////////////////////
		case IloAlgorithm::Feasible:
		{
			UpdateCorePoint(x_tilde);
			AddBendersCuts(x_tilde, theta_tilde);
			if (obj_hat > LowerBound) {
				cout << "Lower bound improved: " << obj_hat << endl;
				UpdateSolution(x_hat, obj_hat);
			}
			if (obj_hat > obj_tilde + EPS) {
				AddTabooCut(model, x_tilde);
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				cout << "Found solution of value: " << obj_hat << endl;
				foundImprovement = true;
				cout << "Found improvement in restricted subproblem." << endl;
			}
			break;
		}
		/////////////////////////////////
		case IloAlgorithm::Optimal:
		{
			UpdateCorePoint(x_tilde);
			AddBendersCuts(x_tilde, theta_tilde);
			AddTrustCutsOuter(x_tilde, 3.0);
			if (obj_hat > LowerBound) {
				cout << "Lower bound improved: " << obj_hat << endl;
				UpdateSolution(x_hat, obj_hat);
			}
			if (obj_hat > obj_tilde + EPS) {	
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				cout << "Found solution of value: " << obj_hat << endl;
				foundImprovement = true;
				cout << "Found improvement and optimal solution in restricted subproblem." << endl;
			}
			break;
		}
		////////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			AddTrustCutsOuter(x_tilde, 3.0);
			cout << "Restricted subproblem was infeasible." << endl;
			break;
		}
		default:
		{
			break;
		}
		} //End of switch statement

	} while (foundImprovement);



	if (diversify) {
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloAlgorithm::Status status = SubproblemDiversify(x_tilde, x_hat, obj_hat);
		switch (status)
		{
			////////////////////////////
		case IloAlgorithm::Feasible:
		{
			cout << "Found feasible solution in diversified subproblem of value: " << obj_hat << endl;
			AddTabooCut(model, x_tilde);
			IntegerCutCallback(x_hat, theta_tilde, obj_hat, false, CallbackStatus::Integer);
			break;
		}
		////////////////////////////
		case IloAlgorithm::Optimal:
		{
			cout << "Found optimal solution in diversified subproblem of value: " << obj_hat << endl;
			AddTrustCutsOuter(x_tilde, 5.0);
			IntegerCutCallback(x_hat, theta_tilde, obj_hat, false, CallbackStatus::Integer);
			break;
		}
		////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			cout << "Diversified subproblem was infeasible." << endl;
			AddTrustCutsOuter(x_tilde, 5.0);
			break;
		}
		////////////////////////////
		default:
		{
			AddTabooCut(model, x_tilde);
			break;
		}
		} //End of switch statement
	}
}

IloAlgorithm::Status Model_LocalBranching::SubproblemRestricted(const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{
	//IloModel submodel(MasterModel);

	IloModel submodel(env);
	submodel.add(model);

	//IloModel submodel = MasterModel;

	IloCplex subcplex(submodel);
	subcplex.use(NULL, 0); //Unregister the callback function


	subcplex.setParam(IloCplex::Param::Tune::Display, 0);
	subcplex.setParam(IloCplex::Param::MIP::Display, 0);
	subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
	subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
	subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
	subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
	subcplex.setParam(IloCplex::Param::Network::Display, 0);
	subcplex.setParam(IloCplex::Param::Conflict::Display, 0);

	subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 3);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 3);

	subcplex.setParam(IloCplex::Param::TimeLimit, 60);

	AddTrustCutsInner(submodel, x_tilde, 2.0);
	//submodel.add(theta_obj >= LowerBound);
	//Two-opt cut generation
	for (int t = 0; t < T; t++) {
		vector<int> always;
		vector<int> active;
		vector<int> inactiveD1;
		vector<int> inactiveD2;
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int k = data.params["station_coord"][j_bar][1];
			if (x_tilde(t, j_bar) > 1 - EPS) {
				int j = data.params["station_coord"][j_bar][0];
				if ((k == ((int)data.params["Mj"][j] - 1)) || (x_tilde(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
				else {
					always.push_back(j_bar);
					submodel.add(x[t][j_bar] == x_tilde(t, j_bar));
				}
			}
			else if ((k == 1) || (x_tilde(t, j_bar - 1) > 1 - EPS)) { inactiveD1.push_back(j_bar); }
			else if ((k == 2) || (x_tilde(t, j_bar - 2) > 1 - EPS)) { inactiveD2.push_back(j_bar); }
			else { submodel.add(x[t][j_bar] == x_tilde(t, j_bar)); }

		}

		//Initialise coverage
		VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
		for (int j_bar : always) { I_tilde += data.a[t].col(j_bar); }
		for (int j_bar : active) { I_tilde += data.a[t].col(j_bar); }

		//Cut for greedy solution + distance 1
		{
			IloExpr lhs(env);
			lhs -= theta[t];
			double covered = (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
			}
			submodel.add(lhs >= -covered);
		}

		//Cuts for adding two outlets (redundant for greedy solution)
		for (int j_bar : inactiveD1) {

			IloExpr lhs(env);
			lhs -= theta[t];

			VectorXd mod_I_tilde = I_tilde + data.a[t].col(j_bar);
			double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
			for (int j_hat : inactiveD2) {
				lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
			}
			submodel.add(lhs >= -covered);
		}

		//Check coverage when swapping out active stations with different ones
		for (int j_bar : active) {

			IloExpr lhs(env);
			lhs -= theta[t];

			VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
			double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
			for (int j_hat : inactiveD1) {
				lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
			}
			submodel.add(lhs >= -covered);
		}
	}

	bool solved = subcplex.solve();
	if (solved) {
		obj = subcplex.getObjValue();
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				sol(t, j_bar) = subcplex.getValue(x[t][j_bar]);
			}
		}
	}
	SolveTime += subcplex.getTime();
	IloAlgorithm::Status status = subcplex.getStatus();
	submodel.end();
	subcplex.end();
	return status;
}

IloAlgorithm::Status Model_LocalBranching::SubproblemDiversify(const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{

	//IloModel submodel(MasterModel);

	IloModel submodel(env);
	submodel.add(model);

	//IloModel submodel = MasterModel;

	IloCplex subcplex(submodel);

	subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
	subcplex.setParam(IloCplex::Param::MIP::Display, 0);
	subcplex.setParam(IloCplex::Param::Tune::Display, 0);
	subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
	subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
	subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
	subcplex.setParam(IloCplex::Param::Network::Display, 0);
	subcplex.setParam(IloCplex::Param::Conflict::Display, 0);

	subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 0);
	subcplex.setParam(IloCplex::Param::Emphasis::MIP, 1);

	subcplex.setParam(IloCplex::Param::TimeLimit, 60);

	//submodel.add(theta_obj >= LowerBound);

	Callback_Improved cb(data, x, theta);
	CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
		| IloCplex::Callback::Context::Id::Relaxation;
	subcplex.use(&cb, contextmask);
	AddTrustCutsInner(submodel, x_tilde, 4.0);
	AddTabooCut(submodel, x_tilde);

	bool solved = subcplex.solve();
	if (solved) {
		obj = subcplex.getObjValue();
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				sol(t, j_bar) = subcplex.getValue(x[t][j_bar]);
			}
		}
	}
	SolveTime += subcplex.getTime();
	IloAlgorithm::Status status = subcplex.getStatus();
	submodel.end();
	subcplex.end();
	return status;
}

double Model_LocalBranching::CalculateObjective(const ArrayXXd& x_tilde)
{
	double covered = 0.0;
	for (int t = 0; t < T; t++) {
		covered += data.Precovered[t];
		VectorXd I_tilde = VectorXd::Zero(P[t]);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) {
				I_tilde += data.a[t].col(j_bar);
				covered += data.Ps[t](j_bar);
			}
		}
		covered += (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
	}

	return covered;
}

vector<VectorXd> Model_LocalBranching::CalculateItilde(const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { I_tilde[t] += data.a[t].col(j_bar); }
		}
	}
	return I_tilde;
}

void Model_LocalBranching::UpdateSolution(const ArrayXXd& x_tilde, double obj)
{
	BestSolution = x_tilde;
	LowerBound = obj;
	thresholdObjective = obj + 1e-4 * (obj + 1e-10);
	UpdatedBound = true;
}

double Model_LocalBranching::GetUpperBound(IloCplex& cplex)
{
	double UB = 0.0;

	IloModel submodel(env);
	submodel.add(model);
	IloCplex subcplex(submodel);

	//vector<IloConversion> conv_x;
	for (int t = 0; t < T; t++) {
		IloConversion conv(env, x[t], IloNumVar::Type::Float);
		submodel.add(conv);
		//conv_x.push_back(conv);
	}

	subcplex.setParam(IloCplex::Param::Emphasis::MIP, 2);
	subcplex.setParam(IloCplex::Param::TimeLimit, 60);
	subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
	subcplex.setParam(IloCplex::Param::Read::WarningLimit, 0);
	subcplex.setParam(IloCplex::Param::MIP::Display, 0);
	subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
	subcplex.use(NULL, 0);
	subcplex.solve();
	UB = subcplex.getBestObjValue();

	subcplex.end();
	submodel.end();

	//vector<IloConversion> conv_x;
	//for (int t = 0; t < T; t++) {
	//	IloConversion conv(env, x[t], IloNumVar::Type::Float);
	//	model.add(conv);
	//	conv_x.push_back(conv);
	//}
	//cplex.use(NULL, 0);
	//cplex.setParam(IloCplex::Param::Emphasis::MIP, 2);
	//cplex.solve();
	//UB = cplex.getBestObjValue();
	//for (auto conv:conv_x) {
	//	model.remove(conv);
	//}
	return UB;
}
