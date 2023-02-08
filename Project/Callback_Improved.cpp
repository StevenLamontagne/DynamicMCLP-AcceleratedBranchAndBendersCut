#include "Callback_Improved.h"


void Callback_Improved::invoke(const IloCplex::Callback::Context& context)
{
	try {
		//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloEnv env = context.getEnv();
		Num2D x_tilde(env, T);
		for (int t = 0; t < T; t++) { x_tilde[t] = IloNumArray(env, M_bar); }

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			double obj_new = context.getCandidateObjective();
			if (obj_new == obj_val) { obj_counter += 1; }
			else { obj_val = obj_new; obj_counter = 0; }
			if (obj_counter >= 25) { context.abort(); }

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				context.getCandidatePoint(x[t], x_tilde[t]);
			}

			LazyCutCallback(context, x_tilde);

			stats["nLazyCuts"] += 1;
			//cout << "Lazy cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;;
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			double obj_new = context.getRelaxationObjective();
			if (obj_new == obj_val) { obj_counter += 1; }
			else { obj_val = obj_new; obj_counter = 0; }
			if (obj_counter >= 25) { context.abort(); }

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				context.getRelaxationPoint(x[t], x_tilde[t]);
			}

			UserCutCallback(context, x_tilde);
			stats["nUserCuts"] += 1;
			//cout << "User cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;

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



void Callback_Improved::LazyCutCallback(const IloCplex::Callback::Context& context, const Num2D& x_tilde)
{
	//IloEnv env = context.getEnv();

	vector<VectorXd> I_tilde;
	for (int t = 0;t<T;t++){ 
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0)); 
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde[t][j_bar] > 1- EPS){I_tilde[t] += data.a[t].col(j_bar);}
		}
	}
	
	AddCuts(context, x_tilde, I_tilde);
	I_tilde.end();
}

void Callback_Improved::UserCutCallback(const IloCplex::Callback::Context& context, const Num2D& x_tilde)
{
	//IloEnv env = context.getEnv();

	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde[t][j_bar] > EPS) { I_tilde[t] += x_tilde[t][j_bar] * data.a[t].col(j_bar); }
		}
	}


	AddCuts(context, x_tilde, I_tilde);
	I_tilde.end();
}

void Callback_Improved::AddCuts(const IloCplex::Callback::Context& context, const Num2D& x_tilde, const vector<VectorXd>& I_tilde)
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
				context.rejectCandidate(lhs >= 0).end();
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -EPS) {
				context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
				//context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			}
			lhs.end();
			break;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	//case multicuts::Multi1PO1:
	//	///////////////////////////////////////////////////////////////////////////////
	//{
	//	for (int t = 0; t < T; t++) {
	//		for (int j = 0; j < M; j++) {
	//			for (int k = 0; k < Mj[j]; k++) {
	//				core_point[t][j][k] = (core_point[t][j][k] + x_tilde[t][j][k]) / 2;
	//			}
	//		}
	//	}

	//	for (int t = 0; t < T; t++) {
	//		IloExpr lhs(env);
	//		lhs -= theta[t][0];
	//		IloNum covered = 0;
	//		bool PO = false;
	//		for (int i = 0; i < N; i++) {
	//			IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
	//			for (int r = 0; r < R[i]; r++) {

	//				if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
	//				else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
	//				else if (I_tilde[t][i][r] > 1) { covered += weight; }
	//				else if (I_tilde[t][i][r] < 1) {
	//					for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
	//						int j = cover_triplet.first;
	//						int k0 = cover_triplet.second;
	//						lhs += weight * x[t][j][k0];
	//					}
	//				}
	//				else {
	//					PO = true;
	//					double I_c = 0;
	//					for (auto cover_triplet : data.cover_triplet[t][i][r]) {
	//						int j = cover_triplet.first;
	//						int k0 = cover_triplet.second;
	//						if (x_tilde[t][j][k0] > EPS) { I_c += x_tilde[t][j][k0]; }
	//					}
	//					if (I_c > 1) { covered += weight; }
	//					else if (I_c < 1) {
	//						for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
	//							int j = cover_triplet.first;
	//							int k0 = cover_triplet.second;
	//							lhs += weight * x[t][j][k0];
	//						}
	//					}
	//					else if (data.P[t][i][r] == triplet::Single) {
	//						int j = data.cover_triplet[t][i][r][0].first;
	//						int k0 = data.cover_triplet[t][i][r][0].second;
	//						lhs += weight * x[t][j][k0];
	//					}
	//					else { covered += weight; }
	//				}

	//			}
	//		}
	//		lhs += covered;
	//		switch (context.getId()) {
	//		case (IloCplex::Callback::Context::Id::Candidate):
	//			if (context.getCandidateValue(lhs) < -EPS) {
	//				if (PO) { stats["nParetoCuts"] += 1; }
	//				context.rejectCandidate(lhs >= 0).end();
	//			}
	//			lhs.end();
	//			break;
	//		case (IloCplex::Callback::Context::Id::Relaxation):
	//			if (context.getRelaxationValue(lhs) < -EPS) {
	//				if (PO) { stats["nParetoCuts"] += 1; }
	//				//context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
	//				context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
	//			}
	//			lhs.end();
	//			break;
	//		}
	//	}

	//	break; //End of switch statement
	//}

}


