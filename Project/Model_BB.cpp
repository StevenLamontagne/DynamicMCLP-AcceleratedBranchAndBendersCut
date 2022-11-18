#include "Model_BB.h"
#define EPS 1.0e-6


ILOSTLBEGIN
typedef IloArray<IloNumArray>    Float2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;

void Model_BB::Solve()
{
	cout << "Creating environment and constants \n";

	time_t start;
	time(&start);
	IloEnv env;

	try {
	//Define constants for easier reading and writing
	int & T = data.T;
	int & M = data.M;
	int & N = data.N;
	vector<int> & Mj = data.Mj;
	vector<int>& R = data.R;


	//Create model and set parameters
	IloModel model(env);
	IloCplex cplex(model);

	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::TimeLimit, 7200);
	cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
	cplex.setParam(IloCplex::Param::WorkMem, 100000);

	//Create variables
	BoolVar3D x(env, T);
	for (int t = 0; t < T; t++) {
		x[t] = BoolVar2D(env,M);
		for (int j = 0; j < M; j++) {
			x[t][j] = IloBoolVarArray(env, Mj[j]);
		}
	}
	
	BoolVar2D y(env, T);
	for (int t = 0; t < T; t++) {
		y[t] = IloBoolVarArray(env, M);
	}

	BoolVar3D w(env, T);
	for (int t = 0; t < T; t++) {
		w[t] = BoolVar2D(env, N);
		for (int i = 0; i < N; i++) {
			w[t][i] = IloBoolVarArray(env, R[i]);
		}
	}


	//Constraints
	////Budget, year 0
	IloExpr budget0(env);
	for (int j = 0; j < M; j++) {
		IloExpr n_outlets(env);
		for (int k = 1; k < Mj[j]; k++) {
			n_outlets += k * x[0][j][k];
		}
		budget0 += (data.params["c"][0][j]) * (n_outlets - data.params["x0"][j]);
		budget0 += ((int)data.params["f"][0][j]) * (y[0][j] - (int)data.params["y0"][j]);
		n_outlets.end();
	}
	model.add(budget0 <= (int)data.params["B"][0]);
	budget0.end();

	////Budget, year 1+
	for (int t = 1; t < T; t++) {
		IloExpr budget(env);
		for (int j = 0; j < M; j++) {
			IloExpr n_outlets(env);
			for (int k = 1; k < Mj[j]; k++) {
				n_outlets += k * (x[t][j][k] - x[t-1][j][k]);
			}
			budget += ((int)data.params["c"][t][j]) * (n_outlets);
			
			budget += ((int)data.params["f"][t][j]) * (y[t][j] - y[t-1][j]);
			n_outlets.end();
		}
		model.add(budget <= (int)data.params["B"][t]);
		budget.end();
	}

	////Can't remove outlets, year 0
	for (int j = 0; j < M; j++) {
		IloExpr KeepX(env);
		int mj = data.params["Mj"][j];
		for (int k = 1; k < mj; k++) {
			KeepX += k * x[0][j][k];
		}
		model.add(KeepX >= (int)data.params["x0"][j]);
		KeepX.end();
	}

	////Can't remove outlets, year 1+
	for (int t = 1; t < T; t++) {
		for (int j = 0; j < M; j++) {
			IloExpr KeepXt1(env);
			IloExpr KeepXt2(env);
			for (int k = 1; k < Mj[j]; k++) {
				KeepXt1 += k * x[t][j][k];
				KeepXt2 += k * x[t-1][j][k];
			}
			model.add(KeepXt1 >= KeepXt2);
			KeepXt1.end();
			KeepXt2.end();
		}
	}

	////Can't close stations, year 0+
	for (int j = 0; j < M; j++) {
		model.add(y[0][j] >= (int)data.params["y0"][j]);
		for (int t = 1; t < T; t++) {
			model.add(y[t][j] >= y[t - 1][j]);
		}
	}

	////Pay one-time cost
	for (int t = 0; t < T; t++) {
		for (int j = 0; j < M; j++) {
			IloExpr open(env);
			for (int k = 1; k < Mj[j]; k++) {
				open += x[t][j][k];
			}
			model.add(open == y[t][j]);
		}
	}

	////Covering 
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < N; i++) {
			for (int r = 0; r < R[i]; r++) {
				switch (data.P[t][i][r])
				{
				case triplet::Uncoverable:
					model.add(w[t][i][r] == 0);
					break;
				case triplet::Precovered:
					model.add(w[t][i][r] == 1);
					break;
				default:
					IloExpr covering(env);
					for (pair<int, int> cover_triplet:data.cover_triplet[t][i][r]){
						int j = cover_triplet.first;
						int k0 = cover_triplet.second;
						for (int k = k0; k < Mj[j]; k++) {
							covering += x[t][j][k];
						}
					}
					model.add(covering >= w[t][i][r]);
					covering.end();
					break;
				}
			}
		}
	}


	//Objective
	IloExpr obj(env);
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < N; i++) {
			obj += ((double)data.params["Ni"][t][i] / (double)data.R[i]) * IloSum(w[t][i]);
		}
	}
	cout << "Model creation time: " << time(NULL) - start << " seconds" << endl;
	model.add(IloMaximize(env, obj));
	obj.end();


	//Solve and get results
	IloBool solved = cplex.solve();
	if (solved) {
		cplex.out() << "Solution status:" << cplex.getCplexStatus() << endl;
		cplex.out() << "Optimal value:" << cplex.getObjValue() << endl;
		for (int t = 0; t < T; t++) {
			vector<int> Solution1;
			for (int j = 0; j < M; j++) {
				int total = 0;
				for (int k = 1; k < Mj[j]; k++) {
					if (cplex.getValue(x[t][j][k]) > 1 - EPS) {
						total = k;
						break;
					}
				}
				Solution1.push_back(total);
			}
			Solution.push_back(Solution1);
		}
		ObjectiveValue = cplex.getObjValue();
		SolveTime = cplex.getTime();
		OptimalityGap = cplex.getMIPRelativeGap();
	}
	else {
		ObjectiveValue = -1;
		SolveTime = -1;
		OptimalityGap = -1;
	}
	}
	catch (IloException & e) {
		cout << "Exception: " << e << endl;
		ObjectiveValue = -1;
		SolveTime = -1;
		OptimalityGap = -1;
	}
	env.end();
}

