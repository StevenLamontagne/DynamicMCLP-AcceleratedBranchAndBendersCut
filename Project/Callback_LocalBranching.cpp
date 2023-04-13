#include "Callback_LocalBranching.h"


//When it's a candidate solution: Run the local branching method
//When it's a fractional solution: Check if there's an integer solution within distance 2 and,
//if so, run the local branching method on it
void Callback_LocalBranching::invoke(const IloCplex::Callback::Context& context)
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
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t,j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Update the core point and the coverage of the core point
void Callback_LocalBranching::UpdateCorePoint(ArrayXXd x_new)
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
We first try and find the best solution within distance 2 of a candidate solution x_tilde using a specialised method. Depending on optimality, we cut the solutions 
around x_tilde. If we can't improve x_tilde, we find a good quality solution near-ish x_tilde to add more diverse Bender's optimality cuts.
*/
void Callback_LocalBranching::IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde, bool diversify)
{
	IloEnv env = context.getEnv();

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
				AddTabooCutMaster(env, context, x_tilde);
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				foundImprovement = true;
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
				AddTrustCutsMaster(env, context, x_tilde, 3.0);
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				foundImprovement = true;
			}
			break;
		}
		////////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			AddTrustCutsMaster(env, context, x_tilde, 3.0);
			break;
		}
		////////////////////////////////
		default:
		{
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
			AddTabooCutMaster(env, context, x_tilde);
			IntegerCutCallback(context, x_hat, obj_hat, false);	
			break;
		}
		////////////////////////////
		case IloAlgorithm::Optimal:
		{
			AddTrustCutsMaster(env, context, x_tilde, 5.0);
			IntegerCutCallback(context, x_hat, obj_hat, false);
			break;
		}
		////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			AddTrustCutsMaster(env, context, x_tilde, 5.0);
			break;
		}
		////////////////////////////
		default:
		{
			AddTabooCutMaster(env, context, x_tilde);
			break;
		}
		} //End of switch statement
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
Add either Multi1B1 or Multi1 Pareto-optimal B1 cuts to the model. Since these are added via rejectCandidate (like the trust cuts), it is not
clear how these are handled by CPLEX. Documentation suggest they are all added as lazy cut constraints, but which ones get added to the model
is unclear.
*/
void Callback_LocalBranching::AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	IloEnv env = context.getEnv();
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
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		}

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
				//context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
				context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
			}
			lhs.end();
			break;
		}
		}
	}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
Add the trust cuts of specified distance to the model. Since we are only considered with integer solutions, the distance is typically 1 higher than we prove
(e.g. if we have solved the subproblem up to distance 2, the constraints will impose a distance of at least 3). Since these constraints require the (fixed size)
trust variables, we replace the cuts by taboo/no-good cuts if we have run out of variables. 
*/
void Callback_LocalBranching::AddTrustCutsMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, double distance)
{
	trustCount += 1; //This should be locked to prevent accidental overwriting, but in single-threaded this is fine

	//We can't change the size of variables midway through solving, and so we must resort to a backup method if we fill up the trust count
	if (trustCount >= trust.getSize()) {
		cout << "Trust cut replaced by no good cut due to variable size limit." << endl;
		AddTabooCutMaster(env, context, x_tilde);
		return;
	}


	IloRangeArray cuts(env);
	for (int t = 0; t < T; t++) {
		IloExpr delta(env);
		delta += distance * trust[trustCount][t];
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
			else { delta += x[t][j_bar]; }
		}
		cuts.add(delta >= distance);
	}
	
	
	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
	{
		context.rejectCandidate(cuts);
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		for (auto i = 0; i < cuts.getSize();i++) {
			context.addUserCut(cuts[i], IloCplex::CutManagement::UseCutForce, false);
		}
		break;
	}
	}

	cuts.end();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
Ensures that the model is solving within the restricted distance. Since this is used for both the restricted and diversified subproblems, this is designed to be added to
any model and with any distance.
*/
void Callback_LocalBranching::AddTrustCutsInner(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, double distance)
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




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
Adds a simple taboo/no-good cut to the model, forcing the solution to have at least one binary variable in x be different than x_tilde.
This version is principally designed for the diversified subproblem, to ensure we end up with a different solution.
*/
void Callback_LocalBranching::AddTabooCut(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde)
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
Adds a simple taboo/no-good cut to the model, forcing the solution to have at least one binary variable in x be different than x_tilde.
This is weaker than the trust counts, but is sometimes the best we can do. Additionally, it does not rely on the trust variables, and so
works as a backup should we run out of those
*/
void Callback_LocalBranching::AddTabooCutMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	IloExpr taboo(env);
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { taboo += 1 - x[t][j_bar]; }
			else { taboo += x[t][j_bar]; }
		}
	}

	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
	{
		context.rejectCandidate(taboo >= 1);
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		context.addUserCut(taboo >= 1, IloCplex::CutManagement::UseCutForce, false);
		break;
	}
	}

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
Checks if the solutions x1 and x2 are within distance 2, using our distance metric.
*/
bool Callback_LocalBranching::WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2)
{
	for (int t = 0; t < T; t++) {
		if ((x1 - x2).abs().sum() > 2.0) { return false; }
	}
	return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
Checks if the solutions x_upper and x are within distance 2, using our distance metric. As the standard WithinRestrictedDistance, but slightly faster
since we do not need the absolute value. 
*/
bool Callback_LocalBranching::WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
{
	for (int t = 0; t < T; t++) {
		if ((x_upper - x).sum() > 2.0) { return false; }
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
Calculates the true objective value for a solution, not relying on the Bender's cuts or CPLEX.
*/
double Callback_LocalBranching::CalculateObjective(const ArrayXXd& x_tilde)
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
Returns the appropriate I_tilde for a given x_tilde. Slightly optimised for integer solutions, since that should be the only
ones which we calculate.
*/
vector<VectorXd> Callback_LocalBranching::CalculateItilde(const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { 
				I_tilde[t] += data.a[t].col(j_bar); }
		}

	}
	return I_tilde;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
Finds the best solution withn distance 2 using a specialised method. More specifically, we generate just enough Bender's optimality cuts
to get the true value of all solutions within distance 2, then run CPLEX (without the callback and with all preprocessing power)
to solve the model directly. In principle, it may be possible to have a slightly faster method using explicit enumeration, however 
this ensures that the managerial constraints (precedence and budget) are satisfied, as well as all of our trust cuts.
In general, this manages to solve the subproblem to optimality.
*/
IloAlgorithm::Status Callback_LocalBranching::SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
Tries to find the best solution within distance 4 (but excluding x_tilde) using branch-and-Bender's cuts. In general, we are not able to solve this and just return 
a heuristic solution different than x_tilde.
*/
IloAlgorithm::Status Callback_LocalBranching::SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
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

	Callback_Improved cb(data, x, theta);
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
