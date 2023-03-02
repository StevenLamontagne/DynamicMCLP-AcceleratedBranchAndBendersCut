#include "Callback_Improved.h"


void Callback_Improved::invoke(const IloCplex::Callback::Context& context)
{
	
	try {
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloEnv env = context.getEnv();
		ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			double obj_new = context.getCandidateObjective();
			if (obj_new == obj_val) { duplicate_counter += 1; }
			else { obj_val = obj_new; duplicate_counter += 0;}
			if (duplicate_counter >= 5)
			{ cout << "Looping detected." << endl; context.pruneCurrentNode(); }
			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}

			LazyCutCallback(context, x_tilde);

			//stats["nLazyCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
			//cout << "Lazy cut time: " << time << endl;;
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			
			if (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) != 0) { break; }
			double obj_new = context.getRelaxationObjective();
			if (obj_new == obj_val) { duplicate_counter += 1; }
			else { obj_val = obj_new; duplicate_counter += 0; }
			if (duplicate_counter >= 5)
			{
				cout << "Looping detected." << endl; context.exitCutLoop();
			}

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}


			UserCutCallback(context, x_tilde);
			stats["nUserCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (UserCutTimes.size() < 512) {UserCutTimes.push_back(time);}
			//cout << "User cut time: " <<  time << endl;

			break;
		}
		default:
			return;
		}
	}
	catch (IloException & e) {
		cout << "Exception during callback: " << e << endl;
		for (pair<string, int> it : stats) {
			if (it.second > 0) {
				std::cout << it.first << " : " << it.second << "\n";
			}
		}

		throw e;
	}
}

void Callback_Improved::UpdateCorePoint(ArrayXXd x_new)
{
	core_point = (core_point + x_new) / 2;
	for (int t = 0; t < T; t++) {
		VectorXd cover = VectorXd::Constant(P[t], 0.0);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (core_point(t,j_bar) > EPS) {
				cover += core_point(t, j_bar) * data.a[t].col(j_bar);
			}
		}
		core_coverage[t] = cover;
	}
}



void Callback_Improved::LazyCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	//IloEnv env = context.getEnv();

	vector<VectorXd> I_tilde;
	for (int t = 0;t<T;t++){ 
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0)); 
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1- EPS){I_tilde[t] += data.a[t].col(j_bar);}
		}
	}
	
	AddCuts(context, x_tilde, I_tilde);
}

void Callback_Improved::UserCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	//IloEnv env = context.getEnv();

	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > EPS) { I_tilde[t] += x_tilde(t, j_bar) * data.a[t].col(j_bar); }
		}
	}


	AddCuts(context, x_tilde, I_tilde);
}

void Callback_Improved::AddCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde)
{
	IloEnv env = context.getEnv();

	//Multi-cut (by year), B1
	for (int t = 0; t < T; t++) {
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		lhs += covered + data.Precovered[t];
		VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
		}

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -EPS) {
				//nCutsAdded += 1;
				stats["nLazyCuts"] += 1;
				context.rejectCandidate(lhs >= 0).end();
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -EPS) {
				//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
				context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			}
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


