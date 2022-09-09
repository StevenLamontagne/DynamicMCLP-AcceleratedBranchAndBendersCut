#include "CoverageCallback.h"
#define EPS 1.0e-6

void CoverageCallback::invoke(const IloCplex::Callback::Context& context)
{
	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate): 
		LazyCutCallback(context);
		break;
	case (IloCplex::Callback::Context::Id::Relaxation):
		UserCutCallback(context);
		break;
	}
}


void CoverageCallback::LazyCutCallback(const IloCplex::Callback::Context& context)
{
	time_t start;
	time(&start);

	IloEnv env = context.getEnv();
	//Get current solution value, x_tilde
	IloArray<IloArray<IloNumArray>> x_tilde(env, T);
	for (int t = 0; t < T; t++) {
		x_tilde[t] = IloArray<IloNumArray>(env, M);
		for (int j = 0; j < M; j++) {
			x_tilde[t][j] = IloNumArray(env, Mj[j]);
			context.getCandidatePoint(x[t][j], x_tilde[t][j]);
		}
	}

	//Calculate I_tilde
	Num3D I_tilde(env, T);
	for (int t = 0; t < T; t++) {
		I_tilde[t] = IloArray<IloArray<IloNum>>(env, N);
		for (int i = 0; i < N; i++) {
			I_tilde[t][i] = IloArray<IloNum>(env, R[i]);
			for (int r = 0; r < R[i]; r++) {
				IloInt val = 0;

				for (int j = 0; j < M; j++) {
					for (int k = 1; k < Mj[j]; k++) {
						if ((1- x_tilde[t][j][k] < EPS) && data.a[t][i][r][j][k]) {
							val += 1;
						}
					}
				}
				I_tilde[t][i][r] = val;
			}
		}
	}
	x_tilde.end();

	AddCuts(context, I_tilde);

}

void CoverageCallback::UserCutCallback(const IloCplex::Callback::Context& context)
{
	time_t start;
	time(&start);

	IloEnv env = context.getEnv();
	
	//Get current solution value, x_tilde
	IloArray<IloArray<IloNumArray>> x_tilde(env, T);
	for (int t = 0; t < T; t++) {
		x_tilde[t] = IloArray<IloNumArray>(env, M);
		for (int j = 0; j < M; j++) {
			x_tilde[t][j] = IloNumArray(env, Mj[j]);
			context.getRelaxationPoint(x[t][j], x_tilde[t][j]);
		}
	}

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
						if ((x_tilde[t][j][k] > EPS) && data.a[t][i][r][j][k]) {
							val +=  x_tilde[t][j][k];
						}
					}
				}

				I_tilde[t][i][r] = val;
			}
		}
	}
	x_tilde.end();

	AddCuts(context, I_tilde);

}

void CoverageCallback::AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde)
{
	IloEnv env = context.getEnv();

	//Generate cuts
	IloExpr lhs(env);
	lhs -= theta;
	IloNum covered = 0;
	switch (cut_type)
	{
	//////////////////////////////////////////////////////////////////////////////////////////
		//Single B0 cut
	case cuts::SingleB0:
	{
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else {
						if (I_tilde[t][i][r] >= 1) {
							covered += weight;
						}
						else {
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		break;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B1 cut
	case cuts::SingleB1:
	{
		for (int t = 0; t < T; t++) {
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
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
								lhs += weight * x[t][j][k];

							}
						}
					}
					else {
						if (I_tilde[t][i][r] >= 1) {
							covered += weight;
						}
						else {
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;

		break;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B2 cut
	case cuts::SingleB2:
	{
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		break;
	}
	//Will probably throw an error earlier, but just in case
	default:
		cout << "Unrecognised cut type" << endl;
		throw;
		break;
	}


	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
		if (context.getCandidateValue(lhs) < -EPS) {
			context.rejectCandidate(lhs >= 0);
		}
		lhs.end();
		break;
	case (IloCplex::Callback::Context::Id::Relaxation):
		if (context.getRelaxationValue(lhs) < -EPS) {
			//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
			context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
		}
		lhs.end();
		break;
	}
}
