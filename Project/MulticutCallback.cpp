#include "MulticutCallback.h"
#define EPS 1.0e-6

void MulticutCallback::invoke(const IloCplex::Callback::Context& context)
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


void MulticutCallback::LazyCutCallback(const IloCplex::Callback::Context& context)
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
						if ((1 - x_tilde[t][j][k] < EPS) && data.a[t][i][r][j][k]) {
							val += 1;
						}
					}
				}
				I_tilde[t][i][r] = val;
			}
		}
	}

	AddCuts(context, I_tilde);
	//context.postHeuristicSolution();
	x_tilde.end();
	I_tilde.end();
}

void MulticutCallback::UserCutCallback(const IloCplex::Callback::Context& context)
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
							val += x_tilde[t][j][k];
						}
					}
				}

				I_tilde[t][i][r] = val;
			}
		}
	}
	

	AddCuts(context, I_tilde);

	x_tilde.end();
	I_tilde.end();

}

void MulticutCallback::AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde)
{
	IloEnv env = context.getEnv();

	switch (cut_type)
	{
	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B1 cut
	case multicuts::SingleB1:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
		break; //End of switch statement
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B2 cut
	case multicuts::SingleB2:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
		break; //End of switch statement
	}
//////////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by year), B1
	case multicuts::Multi1B1:{
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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
			lhs += covered;
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

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by year), B2
	case multicuts::Multi1B2: {
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
			lhs += covered;
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

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by user class), B1
	case multicuts::Multi2B1: {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < T; t++) {
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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
			lhs += covered;
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

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by user class), B2
	case multicuts::Multi2B2: {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < T; t++) {
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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
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
			lhs += covered;
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

		break; //End of switch statement
	}
//////////////////////////////////////////////////////////////////////////////////////////
//Multi-cut (by year and user class), B1
	case multicuts::Multi3B1: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
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
							for (int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
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
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
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
			}
		break; //End of switch statement
	}
////////////////////////////////////////////////////////////////////////////////////
///////Will probably throw an error earlier, but just in case////////////////////////////
	default:
		cout << "Unrecognised cut type" << endl;
		throw;
		break;
	}
/////////////////////////////////////////////////////////////////////////
}
