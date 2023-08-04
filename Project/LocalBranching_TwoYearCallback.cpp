#include "LocalBranching_TwoYearCallback.h"

void LocalBranching_TwoYearCallback::IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde, bool diversify)
{
	IloEnv env = context.getEnv();

	double k_tilde = -1.0;

	//Initially, CPLEX does not have the correct objective for solutions (as the associated Bender's cut is missing). 
	//So, if we receive the solution directly from CPLEX, we calculate the correct objective.
	if (obj_tilde == -1) { obj_tilde = CalculateObjective(x_tilde); }

	ArrayXXd x_hat = ArrayXXd::Zero(T, M_bar);
	double obj_hat = -1.0;

	bool foundImprovement = false;
	do
	{
		foundImprovement = false;
		IloAlgorithm::Status status = SubproblemRestricted(env, x_tilde, x_hat, obj_hat);
		switch (status)
		{
			////////////////////////////////
		case IloAlgorithm::Feasible:
		{
			UpdateCorePoint(x_tilde);
			AddBendersCuts(context, x_tilde);
			if (obj_hat > LowerBound) {
				UpdateSolution(x_hat, obj_hat);
			}

			if (obj_hat > obj_tilde + EPS) {
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				foundImprovement = true;
			}
			else {
				k_tilde = max(k_tilde, 1.0);
			}
			
			break;
		}
		/////////////////////////////////
		case IloAlgorithm::Optimal:
		{
			UpdateCorePoint(x_tilde);
			AddBendersCuts(context, x_tilde);
			if (obj_hat > LowerBound) {
				UpdateSolution(x_hat, obj_hat);
			}

			if (obj_hat > obj_tilde + EPS) {
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				foundImprovement = true;
			}
			else {
				k_tilde = max(k_tilde, 3.0);
			}
			
			break;
		}
		////////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			k_tilde = max(k_tilde, 3.0);
			break;
		}
		////////////////////////////////
		default:
		{
			k_tilde = max(k_tilde, 1.0);
			break;
		}
		} //End of switch statement

	} while (foundImprovement);



	if (diversify) {
		IloAlgorithm::Status status = SubproblemDiversify(env, x_tilde, x_hat, obj_hat);
		switch (status)
		{
			////////////////////////////
		case IloAlgorithm::Feasible:
		{
			if (obj_hat > LowerBound) {
				UpdateSolution(x_hat, obj_hat);
			}
			if (obj_hat > obj_tilde + EPS) {
				IntegerCutCallback(context, x_hat, obj_hat, false);
			}
			break;
		}
		////////////////////////////
		case IloAlgorithm::Optimal:
		{
			if (obj_hat > LowerBound) {
				UpdateSolution(x_hat, obj_hat);
			}

			if (obj_hat > obj_tilde + EPS) {
				IntegerCutCallback(context, x_hat, obj_hat, false);
			}
			else {
				k_tilde = max(k_tilde, 5.0);
			}
			break;
		}
		////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			k_tilde = max(k_tilde, 5.0);
			break;
		}
		////////////////////////////
		default:
		{
			break;
		}
		} //End of switch statement
	}

	//To branch whenever a new solution is found, comment out the first block and
	//uncomment the second block
	//////////////////////////////////////////////////////////////

	//Branch only when the incumbent is updated
	if (obj_tilde > BranchingObjective + EPS) {
		BranchingSolution = x_tilde;
		BranchingObjective = obj_tilde;
		BranchingDistance = k_tilde;
		BranchingFlag = true;
	}

	//////////////////////////////////////////////////////////////

	/////Branch anytime a new solution is found
	//if ((!BranchingFlag) || (obj_tilde > BranchingObjective + EPS)){
	//	BranchingSolution = x_tilde;
	//	BranchingObjective = obj_tilde;
	//	BranchingDistance = k_tilde;
	//	BranchingFlag = true;
	//}

	//////////////////////////////////////////////////////////////
	//end of blocks for branching
}

void LocalBranching_TwoYearCallback::AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	IloEnv env = context.getEnv();
	vector<VectorXd> I_tilde = CalculateItilde(x_tilde);

	
	for (int t = 0; t < T; t++) {

		//To use non-Pareto-optimal B1 cuts, uncomment the first block and comment out the second block
		//////////////////////////////////////////////////////////////

		////Multi-cut (by year), B1
		//	IloExpr lhs(env);
		//	lhs -= theta[t];
		//	double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		//	VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
		//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		//	}

		//////////////////////////////////////////////////////////////

		//Multi-cut (by year), Pareto-optimal B1
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		}

		//////////////////////////////////////////////////////////////
		//end of blocks for cut type

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			if (context.getCandidateValue(lhs) < -covered - EPS) {
				context.rejectCandidate(lhs >= -covered);
			}
			lhs.end();
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			if (context.getRelaxationValue(lhs) < -covered - EPS) {
				context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
			}
			lhs.end();
			break;
		}
		}
	}
}

void LocalBranching_TwoYearCallback::AddTrustCutsInner(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, double distance)
{
	for (int t = 0; t < T; t++) {
		IloExpr delta(env);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
			else { delta += x[t][j_bar]; }
		}
		model.add(delta <= distance);
		delta.end();
	}
}

void LocalBranching_TwoYearCallback::AddTabooCut(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde)
{
	IloExpr taboo(env);
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { taboo += 1 - x[t][j_bar]; }
			else { taboo += x[t][j_bar]; }
		}
	}
	model.add(taboo >= 1);
}

bool LocalBranching_TwoYearCallback::WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2)
{
	for (int t = 0; t < T; t++) {
		if ((x1 - x2).abs().sum() > 2.0) { return false; }
	}
	return true;
}

bool LocalBranching_TwoYearCallback::WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
{
	for (int t = 0; t < T; t++) {
		if ((x_upper - x).sum() > 2.0) { return false; }
	}
	return true;
}

double LocalBranching_TwoYearCallback::CalculateObjective(const ArrayXXd& x_tilde)
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

vector<VectorXd> LocalBranching_TwoYearCallback::CalculateItilde(const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) {
				I_tilde[t] += data.a[t].col(j_bar);
			}
		}

	}
	return I_tilde;
}

void LocalBranching_TwoYearCallback::UpdateCorePoint(ArrayXXd x_new)
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

IloAlgorithm::Status LocalBranching_TwoYearCallback::SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{
	IloModel submodel(env);
	submodel.add(MasterModel);


	IloCplex subcplex(submodel);
	subcplex.use(NULL, 0); //Unregister the callback function, since we use a custom method for distance 2

	subcplex.setParam(IloCplex::Param::Threads, 1);
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

	AddTrustCutsInner(env, submodel, x_tilde, 2.0);

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
	IloAlgorithm::Status status = subcplex.getStatus();
	subcplex.end();
	submodel.end();
	stats["nRestricted"] += 1;

	return status;
}

IloAlgorithm::Status LocalBranching_TwoYearCallback::SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{
	IloModel submodel(env);
	submodel.add(MasterModel);

	IloCplex subcplex(submodel);

	subcplex.setParam(IloCplex::Param::Threads, 1);
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

	MultiCutBenders_Callback cb(data, x, theta);
	CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
		| IloCplex::Callback::Context::Id::Relaxation;
	subcplex.use(&cb, contextmask);
	AddTrustCutsInner(env, submodel, x_tilde, 4.0);
	AddTabooCut(env, submodel, x_tilde);

	bool solved = subcplex.solve();
	if (solved) {
		obj = subcplex.getObjValue();
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				sol(t, j_bar) = subcplex.getValue(x[t][j_bar]);
			}
		}
	}
	IloAlgorithm::Status status = subcplex.getStatus();
	subcplex.end();
	submodel.end();
	stats["nDiversified"] += 1;

	return status;
}

void LocalBranching_TwoYearCallback::invoke(const IloCplex::Callback::Context& context)
{
	try {
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

		switch (context.getId()) {
			////////////////////////////////
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}
			x_tilde = x_tilde.round(); //CPLEX indicates integer solution, so this just eliminates small residuals that could cause problems for the integer callback
			IntegerCutCallback(context, x_tilde);

			stats["nLazyCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }

			break;
		}
		////////////////////////////////
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}

			if (WithinRestrictedDistanceFast(x_tilde, x_tilde.floor())) {
				ArrayXXd x_floor = x_tilde.floor();
				IntegerCutCallback(context, x_floor);
				stats["nUserCuts"] += 1;
			}
			else if (WithinRestrictedDistance(x_tilde, x_tilde.round())) {
				ArrayXXd x_round = x_tilde.round();
				IntegerCutCallback(context, x_round);
				stats["nUserCuts"] += 1;
			}

			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (UserCutTimes.size() < 512) { UserCutTimes.push_back(time); }

			break;
		}
		////////////////////////////////
		case (IloCplex::Callback::Context::Id::Branching):
		{
			if (BranchingFlag) {
				IloEnv env = context.getEnv();

				IloRangeArray cuts_year1(env);
				IloRangeArray cuts_year2(env);

				IloExpr delta_year1(env);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) {
					if (BranchingSolution(0, j_bar) > 1 - EPS) { delta_year1 += 1 - x[0][j_bar]; }
					else { delta_year1 += x[0][j_bar]; }
				}
				cuts_year1.add(delta_year1 >= BranchingDistance);
				cuts_year2.add(delta_year1 <= BranchingDistance - 1);

				IloExpr delta_year2(env);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) {
					if (BranchingSolution(1, j_bar) > 1 - EPS) { delta_year2 += 1 - x[1][j_bar]; }
					else { delta_year2 += x[1][j_bar]; }
				}
				cuts_year2.add(delta_year2 >= BranchingDistance);


				context.makeBranch(cuts_year1, BranchingObjective);
				context.makeBranch(cuts_year2, BranchingObjective);
				stats["nManualBranches"] += 1;
				BranchingFlag = false;
			}
			break;
		}
		////////////////////////////////
		default:
			return;
		}
	}
	catch (IloException & e) {
		cout << "Exception during callback: " << e << endl;
		throw e;
	}

	////Add post-heuristic for best solution found
	////Most values are not correct, but (importantly) the objective is
	////This should allow CPLEX to perform bounds checking
	if (UpdatedSolution) {
		IloEnv env = context.getEnv();

		IloNumVarArray vars(env);
		IloNumArray vals(env);
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				vars.add(x[t][j_bar]);
				vals.add(BestSolution(t, j_bar));
			}
			////CPLEX requires all variables to be included to use the NoCheck option (required to get around trust cuts)
			////So we just set everything to 0
			vars.add(theta[t]);
			vals.add(0);
			for (auto i = 0; i < trust.getSize(); i++) {
				vars.add(trust[i][t]);
				vals.add(0.0);
			}
		}
		vars.add(theta_obj);
		vals.add(LowerBound);

		context.postHeuristicSolution(vars, vals, LowerBound, IloCplex::Callback::Context::SolutionStrategy::NoCheck);
		UpdatedSolution = false;
	}
}
