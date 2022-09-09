#include "MulticutCallback_Strengthened.h"
#define EPS 1.0e-6



void MulticutCallback_Strengthened::invoke(const IloCplex::Callback::Context& context)
{
	//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	IloEnv env = context.getEnv();
	IloArray<IloArray<IloNumArray>> z_tilde(env, T);


	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
	{
		//Get current solution value, z_tilde
		//IloArray<IloArray<IloNumArray>> z_tilde(env, T);
		for (int t = 0; t < T; t++) {
			z_tilde[t] = IloArray<IloNumArray>(env, M);
			for (int j = 0; j < M; j++) {
				z_tilde[t][j] = IloNumArray(env, Mj[j]);
				context.getCandidatePoint(z[t][j], z_tilde[t][j]);
			}
		}

		LazyCutCallback(context, z_tilde);
		//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		//if ((heuristic == useHeuristic::WarmstartAndPostSimple || heuristic == useHeuristic::WarmstartAndPostGreedy) && (context.getIncumbentObjective() != incumbent)) {
		//	incumbent = (double)context.getIncumbentObjective();
		//	AddPostHeuristic(context, z_tilde);
		//}
		//cout << "Lazy cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		//Get current solution value, z_tilde
		//IloArray<IloArray<IloNumArray>> z_tilde(env, T);
		for (int t = 0; t < T; t++) {
			z_tilde[t] = IloArray<IloNumArray>(env, M);
			for (int j = 0; j < M; j++) {
				z_tilde[t][j] = IloNumArray(env, Mj[j]);
				context.getRelaxationPoint(z[t][j], z_tilde[t][j]);
			}
		}


		UserCutCallback(context, z_tilde);
		//if ((heuristic == useHeuristic::WarmstartAndPostSimple || heuristic == useHeuristic::WarmstartAndPostGreedy) && (context.getIncumbentObjective() != incumbent)) {
		//	incumbent = (double)context.getIncumbentObjective();
		//	AddPostHeuristic(context, z_tilde);
		//}
		//cout << "User cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;

		break;
	}
	}
	z_tilde.end();
}


void MulticutCallback_Strengthened::LazyCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& z_tilde)
{
	//time_t start;
	//time(&start);

	IloEnv env = context.getEnv();

	//Calculate I_tilde
	Num3D I_tilde(env, T);
	for (int t = 0; t < T; t++) {
		I_tilde[t] = IloArray<IloArray<IloNum>>(env, N);
		for (int i = 0; i < N; i++) {
			I_tilde[t][i] = IloArray<IloNum>(env, R[i]);
			for (int r = 0; r < R[i]; r++) {
				IloInt val = 0;
				for (auto cover : data.cover[t][i][r]) {
					int j = cover.first;	
					for (int k = 1; k < Mj[j]; k++) {
						if ((data.aBar[t][i][r][j][k] == 1) && (z_tilde[t][j][k] > 1 - EPS)) { val += 1; }
					}
					//int k0 = cover.second;
					//if (z_tilde[t][j][k0] > 1 - EPS) { val += 1; }

					//for (int k = k0; k < Mj[j]; k++) {
					//	if (z_tilde[t][j][k] > 1 - EPS) { val += 1; }
					//}
				}
				I_tilde[t][i][r] = val;
			}
		}
	}


	AddCuts(context, I_tilde);
	I_tilde.end();
}

void MulticutCallback_Strengthened::UserCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& z_tilde)
{
	//time_t start;
	//time(&start);

	IloEnv env = context.getEnv();


	//Calculate I_tilde
	Num3D I_tilde(env, T);
	for (int t = 0; t < T; t++) {
		I_tilde[t] = IloArray<IloArray<IloNum>>(env, N);
		for (int i = 0; i < N; i++) {
			I_tilde[t][i] = IloArray<IloNum>(env, R[i]);
			for (int r = 0; r < R[i]; r++) {
				IloNum val = 0;
				for (int j = 0; j < M; j++) {
					for (int k = 1; k < Mj[j]; k++) {
						if ((z_tilde[t][j][k] > EPS) && data.aBar[t][i][r][j][k]) {
							val += z_tilde[t][j][k];
						}
					}
				}

				I_tilde[t][i][r] = val;
			}
		}
	}


	
	AddCuts(context, I_tilde); 

	I_tilde.end();

}




void MulticutCallback_Strengthened::AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde)
{
	IloEnv env = context.getEnv();

	switch (cut_type)
	{
//////////////////////////////////////////////////////////////////////////////////////////
		//Multi-cut (by year), B1
	case multicuts::Multi1B1: {
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k = data.cover[t][i][r][0].second;
							if (data.aBar[t][i][r][j][k] == 1) { lhs += weight * z[t][j][k]; }
						}
					}
					else {
						if (I_tilde[t][i][r] >= 1) {
							covered += weight;
						}
						else {

							for (pair<int, int> cover : data.cover[t][i][r]) {
								int j = cover.first;
								int k = cover.second;
								if (data.aBar[t][i][r][j][k] == 1) { lhs += weight * z[t][j][k]; }
							}
						}
					}
				}
			}
			lhs += covered;
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

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year and user class), B2
	case multicuts::Multi3B2: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							for (pair<int, int> cover : data.cover[t][i][r]) {
								int j = cover.first;
								int k = cover.second;
								if (data.aBar[t][i][r][j][k] == 1) { lhs += weight * z[t][j][k]; }
							}

						}
					}
				}
				lhs += covered;
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
		}
		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							///////Will probably throw an error earlier, but just in case////////////////////////////
	default:
		cout << "Invalid cut type" << endl;
		throw;
		break;
	}
	/////////////////////////////////////////////////////////////////////////
}

