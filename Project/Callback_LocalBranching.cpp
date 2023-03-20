#include "Callback_LocalBranching.h"



void Callback_LocalBranching::invoke(const IloCplex::Callback::Context& context)
{

	try {

		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloEnv env = context.getEnv();
		ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);


		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			double obj_new = context.getCandidateObjective();
			if (obj_new == obj_val) { obj_counter += 1; }
			else { obj_val = obj_new; obj_counter = 0; }
			if (obj_counter >= 5) {
				//context.pruneCurrentNode();
				context.rejectCandidate();
			} 

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}
			x_tilde = x_tilde.round();
			IntegerCutCallback(context, x_tilde);

			stats["nLazyCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
			//cout << "Lazy cut time: " << time << endl;;
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			double obj_new = context.getRelaxationObjective();
			if ((obj_new - BestObjective) / (1e-10 + BestObjective) > 1e-4) { context.pruneCurrentNode();}

			if ((SkipFractional) || (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) != 0)) { break; }

			
			if (obj_new == obj_val) { obj_counter += 1; }
			else { obj_val = obj_new; obj_counter = 0; }
			if (obj_counter >= 5) {
				context.exitCutLoop();
			}

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t,j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}


			FractionalCutCallback(context, x_tilde);
			stats["nUserCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (UserCutTimes.size() < 512) { UserCutTimes.push_back(time); }
			//cout << "User cut time: " <<  time << endl;

			break;
		}
		default:
			return;
		}
	}
	catch (IloException & e) {
		cout << "Exception during callback: " << e << endl;
		throw e;
	}
	//cout << "Returning to master problem." << endl;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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



void Callback_LocalBranching::IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde, bool diversify)
{
	IloEnv env = context.getEnv();
	
	ArrayXXd x_hat = ArrayXXd::Zero(T, M_bar);
	if (obj_tilde == -1) { obj_tilde = CalculateObjective(x_tilde); }
	double obj_hat = -1.0;

	bool foundImprovement = false;
	do
	{
		foundImprovement = false;

		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloAlgorithm::Status status = SubproblemRestricted(env, x_tilde, x_hat, obj_hat);
		//cout << "Restricted subproblem solving time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		//cout << "Subproblem status: " << status << endl;

		switch (status)
		{
		////////////////////////////////
		case IloAlgorithm::Feasible:
		{
			vector<VectorXd> I_tilde = CalculateItilde(x_tilde);
			UpdateCorePoint(x_tilde);
			AddBendersCuts(context, x_tilde, I_tilde);
			if (obj_hat > BestObjective) {
				cout << "Found solution of value: " << obj_hat << endl;
				BestObjective = obj_hat;
				BestSolution = x_hat;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_hat.row(t).abs() << endl;
				//}
			}
			if (obj_hat > obj_tilde + EPS) {
				//AddTabooCut(env, MasterModel, x_tilde);
				AddTabooCutMaster(env, context, x_tilde);
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				cout << "Found solution of value: " << obj_hat << endl;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_hat.row(t).abs() << endl;
				//}
				foundImprovement = true;
				cout << "Found improvement in restricted subproblem." << endl;
			}
			break;
		}
		/////////////////////////////////
		case IloAlgorithm::Optimal:
		{
			vector<VectorXd> I_tilde = CalculateItilde(x_tilde);
			UpdateCorePoint(x_tilde);
			AddBendersCuts(context, x_tilde, I_tilde);
			if (obj_hat > BestObjective) {
				cout << "Found solution of value: " << obj_hat << endl;
				BestObjective = obj_hat;
				BestSolution = x_hat;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_hat.row(t).abs() << endl;
				//}
			}
			if (obj_hat > obj_tilde + EPS) {
				//AddTrustCutsOuter(env, MasterModel, x_tilde, 3.0);
				AddTrustCutsMaster(env, context, x_tilde, 3.0);
				AddTabooCutMaster(env, context, x_tilde);
				x_tilde = x_hat;
				obj_tilde = obj_hat;
				cout << "Found solution of value: " << obj_hat << endl;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_hat.row(t).abs() << endl;
				//}
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
			//AddTrustCutsOuter(env, MasterModel, x_tilde, 3.0);
			AddTrustCutsMaster(env, context, x_tilde, 3.0);
			AddTabooCutMaster(env, context, x_tilde);
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
		IloAlgorithm::Status status = SubproblemDiversify(env, x_tilde, x_hat, obj_hat);
		//cout << "Diversified subproblem solving time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		//cout << "Subproblem status: " << status << endl;
		switch (status)
		{
		////////////////////////////
		case IloAlgorithm::Feasible:
		{
			cout << "Found feasible solution in diversified subproblem." << endl;
			//AddTabooCut(env, MasterModel, x_tilde);
			AddTabooCutMaster(env, context, x_tilde);
			IntegerCutCallback(context, x_hat, obj_hat, false);	
			break;
		}
		////////////////////////////
		case IloAlgorithm::Optimal:
		{
			cout << "Found optimal solution in diversified subproblem." << endl;
			//AddTrustCutsOuter(env, MasterModel, x_tilde,  5.0);
			AddTrustCutsMaster(env, context, x_tilde, 5.0);
			IntegerCutCallback(context, x_hat, obj_hat, false);
			break;
		}
		////////////////////////////
		case IloAlgorithm::Infeasible:
		case IloAlgorithm::Unbounded:
		case IloAlgorithm::InfeasibleOrUnbounded:
		{
			cout << "Diversified subproblem was infeasible." << endl;
			//AddTrustCutsOuter(env, MasterModel, x_tilde, 5.0);
			AddTrustCutsMaster(env, context, x_tilde, 5.0);
			break;
		}
		////////////////////////////
		default:
		{
			//AddTabooCut(env, MasterModel, x_tilde);
			AddTabooCutMaster(env, context, x_tilde);
			break;
		}
		} //End of switch statement
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Callback_LocalBranching::FractionalCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde)
{
	//IloEnv env = context.getEnv();
	if (WithinRestrictedDistanceFast(x_tilde, x_tilde.floor())) {
		ArrayXXd x_floor = x_tilde.floor();
		//IntegerCutCallback(context, x_floor);
	}
	else if (WithinRestrictedDistance(x_tilde, x_tilde.round())) {
		ArrayXXd x_round = x_round.floor();
		//IntegerCutCallback(context, x_round);
	}
	else{
		vector<VectorXd> I_tilde;
		for (int t = 0; t < T; t++) {
			I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				if (x_tilde(t, j_bar) > EPS) { I_tilde[t] += x_tilde(t, j_bar) * data.a[t].col(j_bar); }
			}
		}
		AddBendersCuts(context, x_tilde, I_tilde);
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Callback_LocalBranching::AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde)
{
	IloEnv env = context.getEnv();

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
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		}

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -covered - EPS) {
				context.rejectCandidate(lhs >= -covered);
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -covered - EPS) {
				//context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
				context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
			}
			lhs.end();
			break;
		}
	}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Callback_LocalBranching::AddTrustCutsOuter(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, double distance)
{
	trustCount += 1;

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

	model.add(IloSum(trust[trustCount]) <= T - 1);



}

void Callback_LocalBranching::AddTrustCutsMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, double distance)
{
	trustCount += 1;
	IloRangeArray cuts(env);
	cuts.add(IloSum(trust[trustCount]) <= T - 1);
	//cout << cuts[0] << endl;
	for (int t = 0; t < T; t++) {
		IloExpr delta(env);
		delta += distance * trust[trustCount][t];
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
			else { delta += x[t][j_bar]; }
		}
		cuts.add(delta >= distance);
		//cout << cuts[t+1] << endl;
		delta.end();
	}
	
	context.rejectCandidate(cuts);

	//context.rejectCandidate(IloSum(trust[trustCount]) <= T - 1);
	//for (int t = 0; t < T; t++) {
	//	IloExpr delta(env);
	//	delta += distance * trust[trustCount][t];
	//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
	//		if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
	//		else { delta += x[t][j_bar]; }
	//	}
	//	context.rejectCandidate(delta >= distance);
	//	delta.end();
	//}

	//cuts.end();

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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

void Callback_LocalBranching::AddTabooCutMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	IloExpr taboo(env);
	for (int t = 0; t < T; t++) {
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { taboo += 1 - x[t][j_bar]; }
			else { taboo += x[t][j_bar]; }
		}
	}
	//cout << taboo << endl;
	context.rejectCandidate(taboo >= 1);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////Calculate the l-infinity norm across all years for each j_bar, then return the sum 
//double Callback_LocalBranching::CalculateDistance(const ArrayXXd& x1, const ArrayXXd& x2)
//{
//	ArrayXd distance = ArrayXd::Zero(M_bar);
//	for (int t = 0; t < T; t++) {
//		distance.max((x1 - x2).abs());
//	}
//	return distance.sum();
//}
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////As the CalculateDistance function, but when we are guaranteed that x_upper >= x 
//double Callback_LocalBranching::CalculateDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
//{
//	ArrayXd distance = ArrayXd::Zero(M_bar);
//	for (int t = 0; t < T; t++) {
//		distance.max(x_upper - x);
//	}
//	return distance.sum();
//}


//Calculate the l-infinity norm across all years for each j_bar, then return the sum 
bool Callback_LocalBranching::WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2)
{
	for (int t = 0; t < T; t++) {
		if ((x1 - x2).abs().sum() > 2.0) { return false; }
	}
	return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//As the CalculateDistance function, but when we are guaranteed that x_upper >= x 
bool Callback_LocalBranching::WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
{
	for (int t = 0; t < T; t++) {
		if ((x_upper - x).sum() > 2.0) { return false; }
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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


vector<VectorXd> Callback_LocalBranching::CalculateItilde(const ArrayXXd& x_tilde)
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



IloAlgorithm::Status Callback_LocalBranching::SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{
	//IloModel submodel(MasterModel);

	IloModel submodel(env);
	submodel.add(MasterModel);

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
				if ((k == ((int) data.params["Mj"][j] - 1)) || (x_tilde(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
				else { 
					always.push_back(j_bar); 
					submodel.add(x[t][j_bar] == x_tilde(t, j_bar));
				}
			}
			else if ((k == 1) || (x_tilde(t, j_bar - 1) > 1 - EPS)) { inactiveD1.push_back(j_bar); }
			else if ((k == 2) || (x_tilde(t, j_bar - 2) > 1 - EPS)) { inactiveD2.push_back(j_bar); }
			else{ submodel.add(x[t][j_bar] == x_tilde(t, j_bar)); }

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
			double covered= (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
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
				sol(t,j_bar) = subcplex.getValue(x[t][j_bar]);
			}
		}
	}

	//IloAlgorithm::Status
	return subcplex.getStatus();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



IloAlgorithm::Status Callback_LocalBranching::SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj)
{

	
	//IloModel submodel(MasterModel);

	IloModel submodel(env);
	submodel.add(MasterModel);

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

	return subcplex.getStatus();
}
