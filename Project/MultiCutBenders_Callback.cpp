#include "MultiCutBenders_Callback.h"


void MultiCutBenders_Callback::invoke(const IloCplex::Callback::Context& context)
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
			else { obj_val = obj_new; duplicate_counter = 0;}
			if (duplicate_counter >= 5) { cout << "Looping detected at integer node." << endl; 
			//context.rejectCandidate(); 
			break;}

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}

			LazyCutCallback(context, x_tilde);
			UpdateCorePoint(x_tilde);

			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			
			if (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) != 0) { break; }
			double obj_new = context.getRelaxationObjective();
			if (obj_new == obj_val) { duplicate_counter += 1; }
			else { obj_val = obj_new; duplicate_counter = 0; }
			if (duplicate_counter >= 5){ cout << "Looping detected at fractional node." << endl; 
			//context.exitCutLoop(); 
			break; }

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}


			UserCutCallback(context, x_tilde);

			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (UserCutTimes.size() < 512) {UserCutTimes.push_back(time);}


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

void MultiCutBenders_Callback::UpdateCorePoint(ArrayXXd x_new)
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



void MultiCutBenders_Callback::LazyCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0;t<T;t++){ 
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0)); 
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1- EPS){I_tilde[t] += data.a[t].col(j_bar);}
		}
	}
	
	AddCuts(context, x_tilde, I_tilde);
}

void MultiCutBenders_Callback::UserCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > EPS) { I_tilde[t] += x_tilde(t, j_bar) * data.a[t].col(j_bar); }
		}
	}


	AddCuts(context, x_tilde, I_tilde);
}

void MultiCutBenders_Callback::AddCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde)
{
	IloEnv env = context.getEnv();

	
	for (int t = 0; t < T; t++) {
		IloExpr lhs(env);
		lhs -= theta[t];

		//Multi-cut (by time period), B1
		//double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		//VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
		//for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//	lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		//}


		//Multi-cut (by time period), B1 + pareto-optimal
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);
		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		}


		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -covered - EPS) {
				context.rejectCandidate(lhs >= -covered);
				stats["nLazyCuts"] += 1;

			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -covered - EPS) {
				context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
				stats["nUserCuts"] += 1;
			}
			lhs.end();
			break;
		}
	}





}


