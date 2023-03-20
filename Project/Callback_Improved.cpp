#include "Callback_Improved.h"


void Callback_Improved::invoke(const IloCplex::Callback::Context& context)
{
	
	try {
		chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		IloEnv env = context.getEnv();
		ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
		ArrayXd theta_tilde = ArrayXd::Constant(T, 0.0);

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{


			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				theta_tilde(t) = context.getCandidatePoint(theta[t]);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
			}
			
			
			double obj_new = context.getCandidateObjective();

//			if (((theta_tilde == IncumbentTheta).all()) && ((x_tilde == IncumbentSolution).all())) {
			bool same_theta = false;
			bool same_x = false;
			if (obj_new == IncumbentObjective) {
				same_theta = true;
			}
			if ( (same_theta) && ((x_tilde == IncumbentSolution).all())) {
				same_x = true;	
				//cout << " x variables have not changed." << endl;
			}
			if (same_theta && same_x) {
				//cout << " Theta and x variables have not changed." << endl;
				duplicate_counter += 1;
			}
			else { 
				IncumbentTheta = theta_tilde; 
				IncumbentSolution = x_tilde; 
				IncumbentObjective = obj_new;
				duplicate_counter = 0;
			}
			//double obj = data.SolutionQuality(x_tilde);
			//if (obj > BestObjective) {
			//	BestObjective = obj;
			//	cout << "Found solution of value: " << obj << endl;
			//}

			if (duplicate_counter >= 5)
			{
				cout << "Looping detected at integer node." << endl;

				//cout << "Objective value: " << obj_new << endl;
				//cout << "Exact value: " << data.SolutionQuality(x_tilde) << endl;
				//cout << "Solution: " << endl;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_tilde.row(t) << endl;
				//}
				//cout << "Theta: " << theta_tilde << endl;

				vector<VectorXd> I_tilde;
				for (int t = 0; t < T; t++) {
					I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
					for (int j_bar = 0; j_bar < M_bar; j_bar++) {
						if (x_tilde(t, j_bar) > 1 - EPS) { 
							I_tilde[t] += data.a[t].col(j_bar); 
						}
					}
				}
				//

				for (int t = 0; t < T; t++) {
					IloExpr lhs(env);
					lhs -= theta[t];
					
					//double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
					//VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();

					double covered= ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);
					VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
					
					for (int j_bar = 0; j_bar < M_bar; j_bar++) {
						lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
					}
					double val = context.getCandidateValue(lhs);
					if (val < -covered- EPS) {
						cout << "Cut value in year " << t << ": " << val << ", violation amount: " << covered + val << endl;
					}
					
				}

				duplicate_counter = 0;
				context.rejectCandidate();
			}
			UpdateCorePoint(x_tilde);
			LazyCutCallback(context, x_tilde);

			stats["nLazyCuts"] += 1;
			double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			
			if ((SkipFractional) || (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) != 0)) { break; }

			//Get current solution value, x_tilde and theta_tilde
			for (int t = 0; t < T; t++) {
				theta_tilde(t) = context.getRelaxationPoint(theta[t]);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}

			double obj_new = context.getRelaxationObjective();
			//if (((theta_tilde == IncumbentTheta).all()) && ((x_tilde == IncumbentSolution).all())) { duplicate_counter += 1; }
			bool same_theta = false;
			bool same_x = false;
			if (obj_new == IncumbentObjective) {
				same_theta = true;
			}
			//if ((theta_tilde == IncumbentTheta).all()) {
			//	same_theta = true;

			//}
			if ((same_theta) && ((x_tilde == IncumbentSolution).all())) {
				same_x = true;
				//cout << " x variables have not changed." << endl;
			}
			if (same_theta && same_x) {
				//cout << " Theta and x variables have not changed." << endl;
				duplicate_counter += 1;
			}
			else {
				IncumbentTheta = theta_tilde;
				IncumbentSolution = x_tilde;
				IncumbentObjective = obj_new;
				duplicate_counter = 0;
			}
			if (duplicate_counter >= 5)
			{
				cout << "Looping detected at fractional node." << endl;
				//cout << "Objective value: " << obj_new << endl;
				//cout << "Exact value: " << data.SolutionQualityContinuous(x_tilde) << endl;

				//cout << "Solution: " << endl;
				//for (int t = 0; t < T; t++) {
				//	cout << "Year " << t << ": " << x_tilde.row(t) << endl;
				//}
				//cout << "Theta: " << theta_tilde << endl;

				SkipFractional = true;
				vector<VectorXd> I_tilde;
				for (int t = 0; t < T; t++) {
					I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
					for (int j_bar = 0; j_bar < M_bar; j_bar++) {
						if (x_tilde(t, j_bar) > EPS) { I_tilde[t] += x_tilde(t, j_bar) * data.a[t].col(j_bar); }
					}
				}


				for (int t = 0; t < T; t++) {
					IloExpr lhs(env);
					lhs -= theta[t];

					//double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
					//VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();

					double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);
					VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();

					for (int j_bar = 0; j_bar < M_bar; j_bar++) {
						lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
					}
					double val = context.getRelaxationValue(lhs);
					context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
					if (val < -covered - EPS) {
						cout << "Cut value in year " << t << ": " << val << ", violation amount: " << covered + val << endl;
					}

				}
				context.exitCutLoop();
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
			if (context.getCandidateValue(lhs) < -covered -EPS) {
				context.rejectCandidate(lhs >= -covered);
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -covered -EPS) {
				context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
				//context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
			}
			lhs.end();
			break;
		}
	}

}


