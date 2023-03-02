#include "Callback_LocalBranching.h"



void Callback_LocalBranching::invoke(const IloCplex::Callback::Context& context)
{

	try {
		double obj_new = context.getCandidateObjective();
		if (obj_new == obj_val) { obj_counter += 1; }
		else { obj_val = obj_new; obj_counter = 0; }

		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloEnv env = context.getEnv();
		ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);


		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			if (obj_counter >= 25) {
				//context.pruneCurrentNode();
				context.exitCutLoop();
			} 

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}
			double obj_tilde = context.getCandidateObjective();
			IntegerCutCallback(context, x_tilde, obj_tilde);

			stats["nLazyCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
			//cout << "Lazy cut time: " << time << endl;;
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			if (obj_counter >= 25) {
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
	double obj_hat = -1.0;

	bool foundImprovement = false;
	do
	{
		foundImprovement = false;
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloAlgorithm::Status status = SubproblemRestricted(env, x_tilde, x_hat, obj_hat);
		cout << "Restricted subproblem solving time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		cout << "Subproblem status: " << status << endl;

		UpdateCorePoint(x_tilde);
		vector<VectorXd> I_tilde = CalculateItilde(x_tilde);
		AddBendersCuts(context, x_tilde, I_tilde);
		int temp = 0;

		if ((status == IloAlgorithm::Status::Infeasible) || (status == IloAlgorithm::Status::InfeasibleOrUnbounded)) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 3.0);
			obj_hat = -1.0;
			cout << "Restricted subproblem was infeasible." << endl;
		}
		else if ((obj_hat > obj_tilde + EPS) && (status == IloAlgorithm::Status::Optimal)) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 3.0);
			x_tilde = x_hat;
			obj_tilde = obj_hat;
			foundImprovement = true;
			cout << "Found improvement and optimal solution in restricted subproblem." << endl;
		}
		else if ((obj_hat > obj_tilde + EPS) && (status == IloAlgorithm::Status::Feasible)) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 2.0);
			x_tilde = x_hat;
			obj_tilde = obj_hat;
			foundImprovement = true;
			cout << "Found improvement and feasible solution in restricted subproblem." << endl;
		}

	} while (foundImprovement);

	if (diversify) {
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloAlgorithm::Status status = SubproblemDiversify(env, x_tilde, x_hat, obj_hat);
		cout << "Diversified subproblem solving time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		cout << "Subproblem status: " << status << endl;
		if ((status == IloAlgorithm::Status::Infeasible) || (status == IloAlgorithm::Status::InfeasibleOrUnbounded)) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 5.0);
			//cout << "Diversified subproblem was infeasible." << endl;
		}
		else if (status == IloAlgorithm::Status::Optimal) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 5.0);
			IntegerCutCallback(context, x_hat, obj_hat, false);
			//cout << "Found optimal solution in diversified subproblem." << endl;
		}
		else if (status == IloAlgorithm::Status::Feasible) {
			AddTrustCuts(env, *MasterModel, x_tilde, -1.0, 3.0);
			IntegerCutCallback(context, x_hat, obj_hat, false);
			//cout << "Found feasible solution in diversified subproblem." << endl;
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Callback_LocalBranching::FractionalCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde)
{
	//IloEnv env = context.getEnv();
	if (CalculateDistanceFast(x_tilde, x_tilde.floor()) <= 2.0) {
		ArrayXXd x_floor = x_tilde.floor();
		double covered = CalculateObjective(x_floor);
		IntegerCutCallback(context, x_floor, covered);
	}
	else if (CalculateDistance(x_tilde, x_tilde.round()) <= 2.0) {
		ArrayXXd x_round = x_round.floor();
		double covered = CalculateObjective(x_round);
		IntegerCutCallback(context, x_round, covered);
	}
	else if (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) == 0) {
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
	
	//Multi-cut (by year), B1
	for (int t = 0; t < T; t++) {
		double obj = 0.0;
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		lhs += covered + data.Precovered[t];
		obj += covered + data.Precovered[t];
		VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
			if (x_tilde(t, j_bar) > EPS) {
				obj += data.Ps[t](j_bar);
			}
		}
		obj += (I_tilde[t].array() < 1).select(I_tilde[t], VectorXd::Zero(P[t])).matrix().cast<double>().dot(data.weights[t]);

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (obj < context.getCandidateObjective()) {
				context.rejectCandidate(lhs >= 0).end();
			}
			//if (context.getCandidateValue(lhs) < -EPS) {
			//	context.rejectCandidate(lhs >= 0).end();
			//}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (obj < context.getRelaxationObjective()) {
				context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			}
			//if (context.getRelaxationValue(lhs) < -EPS) {
			//	context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
			//	//context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			//}
			lhs.end();
			break;
		}
	}



	////Multi-cut (by year), Multi1PO1
	//for (int t = 0; t < T; t++) {
	//	IloExpr lhs(env);
	//	lhs -= theta[t];
	//	double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

	//	lhs += covered + data.Precovered[t];
	//	VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
	//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
	//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
	//	}

	//	switch (context.getId()) {
	//	case (IloCplex::Callback::Context::Id::Candidate):
	//		if (context.getCandidateValue(lhs) < -EPS) {
	//			context.rejectCandidate(lhs >= 0).end();
	//		}
	//		lhs.end();
	//		break;
	//	case (IloCplex::Callback::Context::Id::Relaxation):
	//		if (context.getRelaxationValue(lhs) < -EPS) {
	//			context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
	//			//context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
	//		}
	//		lhs.end();
	//		break;
	//	}
	//}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Callback_LocalBranching::AddTrustCuts(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, const double& minimum, const double& maximum)
{

	IloNumVarArray delta_bar(env, M_bar, 0.0, 1.0);
	BoolVar2D offset(env, M_bar); //Note that the usual order is inverted here (i.e. offset[j_bar][t] instead of offset[t][j_bar])


	//Create variables and constraints for calculating the l-infinity norm across all years for each j_bar
	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		offset[j_bar] = IloBoolVarArray(env, T);
		model.add(IloSum(offset[j_bar]) <= T - 1);
		for (int t = 0; t < T; t++) {
			IloExpr distance;
			if (x_tilde(t, j_bar) > 1 - EPS) { distance = 1 - x[t][j_bar]; }
			else { distance = x[t][j_bar]; }

			model.add(delta_bar[j_bar] >= distance);
			model.add(delta_bar[j_bar] - offset[j_bar][t] <= distance);
			distance.end();
		}
	}

	//Minimum should always be greater than zero. Otherwise, should use the other version (skips adding binary variables)
	model.add(IloSum(delta_bar) >= minimum);

	//Maximum distance
	if (maximum >= 0) {
		model.add(IloSum(delta_bar) <= maximum);
	}


}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Callback_LocalBranching::AddTrustCuts(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, const double& maximum)
{
	IloNumVarArray delta_bar(env, M_bar, 0.0, 1.0);

	//Create constraints for calculating the l-infinity norm across all years for each j_bar
	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		for (int t = 0; t < T; t++) {
			IloExpr distance;
			if (x_tilde(t, j_bar) > 1 - EPS) { distance = 1 - x[t][j_bar]; }
			else { distance = x[t][j_bar]; }

			model.add(delta_bar[j_bar] >= distance);
			distance.end();
		}
	}

	//Maximum should always be greater than zero. Otherwise, should use the other version (needs binary variables)
	model.add(IloSum(delta_bar) <= maximum);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Calculate the l-infinity norm across all years for each j_bar, then return the sum 
double Callback_LocalBranching::CalculateDistance(const ArrayXXd& x1, const ArrayXXd& x2)
{
	ArrayXd distance = ArrayXd::Zero(M_bar);
	for (int t = 0; t < T; t++) {
		distance.max((x1 - x2).abs());
	}
	return distance.sum();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//As the CalculateDistance function, but when we are guaranteed that x_upper >= x 
double Callback_LocalBranching::CalculateDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
{
	ArrayXd distance = ArrayXd::Zero(M_bar);
	for (int t = 0; t < T; t++) {
		distance.max(x_upper - x);
	}
	return distance.sum();
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
	IloModel submodel(*MasterModel);
	IloCplex subcplex(submodel);
	subcplex.setParam(IloCplex::Param::TimeLimit, 60);
	subcplex.setParam(IloCplex::Param::MIP::Display, 0);
	subcplex.setParam(IloCplex::Param::Tune::Display, 0);
	subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
	subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
	subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
	subcplex.setParam(IloCplex::Param::Network::Display, 0);
	subcplex.setParam(IloCplex::Param::Conflict::Display, 0);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 0);

	subcplex.use(NULL, 0); //Unregister the callback function
	AddTrustCuts(env, submodel, x_tilde, 2.0);

	//Two-opt cut generation
	for (int t = 0; t < T; t++) {
		vector<int> always;
		vector<int> active;
		vector<int> inactive;
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int k = data.params["station_coord"][j_bar][1];
			if (x_tilde(t, j_bar) > 1 - EPS) {
				int j = data.params["station_coord"][j_bar][0];
				if ((k == data.params["Mj"][j] - 1) || (x_tilde(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
				else { always.push_back(j_bar); }
			}
			else if ((k <= 2) || (x_tilde(t, j_bar - 1) > 1 - EPS) || ((k >= 3) && (x_tilde(t, j_bar - 2) > 1 - EPS))) { inactive.push_back(j_bar); }
		}
		//for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//	int k = data.params["station_coord"][j_bar][1];
		//	if (x_tilde(t, j_bar) > 1 - EPS) {
		//		int j = data.params["station_coord"][j_bar][0];
		//		if ((k == data.params["Mj"][j] - 1) || (x_tilde(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
		//		else { always.push_back(j_bar); }
		//	}
		//	else if ((k == 1) || (x_tilde(t, j_bar - 1) > 1 - EPS)) { inactive.push_back(j_bar); }
		//}

		//Initialise coverage
		VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
		for (int j_bar : always) { I_tilde += data.a[t].col(j_bar); }
		for (int j_bar : active) { I_tilde += data.a[t].col(j_bar);}

		//Cut for greedy solution
		{ 
			IloExpr lhs(env);
			lhs -= theta[t];
			lhs += data.Precovered[t];

			double covered = (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			lhs += covered + data.Precovered[t];
			VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
			}
			submodel.add(lhs >= 0);
		}

		
		//Initialise covered (does not include active set)
		double covered = data.Precovered[t];
		for (int j_bar : always) {
			covered += data.Ps[t](j_bar);
		}

		//Check coverage when swapping out active stations with different ones
		for (int j_bar : active) {
			//IloExpr lhs(env);
			//lhs -= theta[t];

			//VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
			//double mod_covered = covered;
			//mod_covered += (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			//for (int j_hat : active) {
			//	if (j_hat == j_bar) { continue; }
			//	mod_covered += data.Ps[t](j_hat);
			//}
			//lhs += mod_covered;


			//VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
			//lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
			//for (int j_hat : inactive) {
			//	lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered) + data.Ps[t](j_hat)) * x[t][j_hat];
			//}
			//submodel.add(lhs >= 0);


			IloExpr lhs(env);
			lhs -= theta[t];

			VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
			double mod_covered = covered;
			mod_covered += (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			for (int j_hat : active) {
				lhs += data.Ps[t](j_hat) * x[t][j_hat];
			}
			lhs += mod_covered;


			VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
			for (int j_hat : inactive) {
				lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered) + data.Ps[t](j_hat)) * x[t][j_hat];
			}
			submodel.add(lhs >= 0);
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

	
	IloModel submodel(*MasterModel);
	IloCplex subcplex(submodel);
	subcplex.setParam(IloCplex::Param::TimeLimit, 60);
	//subcplex.setParam(IloCplex::Param::MIP::Display, 0);
	subcplex.setParam(IloCplex::Param::Tune::Display, 0);
	//subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
	subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
	subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
	subcplex.setParam(IloCplex::Param::Network::Display, 0);
	subcplex.setParam(IloCplex::Param::Conflict::Display, 0);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 0);
	subcplex.setParam(IloCplex::Param::Emphasis::MIP, 1);

	Callback_Improved cb(data, x, theta);
	CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
		| IloCplex::Callback::Context::Id::Relaxation;
	subcplex.use(&cb, contextmask);
	AddTrustCuts(env, submodel, x_tilde, 3.0, 4.0);

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
