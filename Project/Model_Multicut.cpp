#include "Model_Multicut.h"
#define EPS 1.0e-6

ILOSTLBEGIN
typedef IloArray<IloNumArray>    Float2D;
typedef IloArray<IloNumVarArray>    FloatVar2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;


void Model_Multicut::Solve(multicuts cut_type, vector<vector<int>>* warmstart)
{
	time_t start;
	time(&start);
	IloEnv env;

	try {
		//Define constants for easier reading and writing
		int& T = data.T;
		int& M = data.M;
		int& N = data.N;
		vector<int>& Mj = data.Mj;
		vector<int>& R = data.R;


		//Create model and set parameters
		IloModel model(env);
		IloCplex cplex(model);

		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200);
		cplex.setParam(IloCplex::Param::Preprocessing::Linear, 0);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);

		//Create variables
		BoolVar3D x(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = BoolVar2D(env, M);
			for (int j = 0; j < M; j++) {
				x[t][j] = IloBoolVarArray(env, Mj[j]);
			}
		}

		BoolVar2D y(env, T);
		for (int t = 0; t < T; t++) {
			y[t] = IloBoolVarArray(env, M);
		}

		FloatVar2D theta(env, T);
		for (int t = 0; t < T; t++) {
			theta[t] = IloNumVarArray(env, N, 0.0, IloInfinity);
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
					n_outlets += k * (x[t][j][k] - x[t - 1][j][k]);
				}
				budget += ((int)data.params["c"][t][j]) * (n_outlets);

				budget += ((int)data.params["f"][t][j]) * (y[t][j] - y[t - 1][j]);
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
					KeepXt2 += k * x[t - 1][j][k];
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

		///Redundant, but prevents crash from unextracted variables
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				model.add(x[t][j][0] >= 0);
			}
		}


		//Cut specific constraints and objective
		/////////////////////////////////////////////////////////////////////////
		switch (cut_type)
		{
			///////////////////////////////////////////////////////////
		case multicuts::SingleB1:
		case multicuts::SingleB2:
		{
			////Upper bound for theta (ensures bounded problem)
			IloNum ub = 0;
			for (int t = 0; t < T; t++) {
				for (int i = 0; i < N; i++) {
					ub += (IloInt)data.params["Ni"][t][i];
				}
			}
			model.add(theta[0][0] <= ub);


			//Objective
			cout << "Adding objective \n";
			IloExpr obj(env);
			obj += theta[0][0];
			model.add(IloMaximize(env, obj));
			obj.end();

			//Warmstart
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);

			if (!warmstart) {
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						startVal.add((int)data.params["y0"][j]);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (data.params["x0"][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
				}
				startVar.add(theta[0][0]);
				startVal.add(0);
			}
			else {
				auto sol = *warmstart;
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						int val = 0;
						if (sol[t][j] > 0) { val = 1; }
						startVal.add(val);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (sol[t][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
				}
				startVar.add(theta[0][0]);
				startVal.add(data.SolutionQuality(sol));
			}

			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			break;
		}
			///////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////
		case multicuts::Multi1B1:
		case multicuts::Multi1B2:
		{
			////Upper bound for theta (ensures bounded problem)

			for (int t = 0; t < T; t++) {
				IloNum ub = 0;
				for (int i = 0; i < N; i++) {
					ub += (IloInt)data.params["Ni"][t][i];
				}
				model.add(theta[t][0] <= ub);
			}



			//Objective
			cout << "Adding objective \n";
			IloExpr obj(env);
			for (int t = 0; t < T; t++) { obj += theta[t][0]; }
			model.add(IloMaximize(env, obj));
			obj.end();

			//Warmstart
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);

			if (!warmstart) {
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						startVal.add((int)data.params["y0"][j]);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (data.params["x0"][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
					startVar.add(theta[t][0]);
					startVal.add(0);
				}
			}
			else {
				auto sol = *warmstart;
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						int val = 0;
						if (sol[t][j] > 0) { val = 1; }
						startVal.add(val);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (sol[t][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
					startVar.add(theta[t][0]);
					IloNum val = 0;
					for (int i = 0; i < N; i++) {
						IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
						for (int r = 0; r < R[i]; r++) {
							if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						}
					}
					startVal.add(val);
				}

			}

			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			break;
		}
			///////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////
		case multicuts::Multi2B1:
		case multicuts::Multi2B2:
		{
			////Upper bound for theta (ensures bounded problem)

			
				
			for (int i = 0; i < N; i++) {
				IloNum ub = 0;
				for (int t = 0; t < T; t++) {
					ub += (IloInt)data.params["Ni"][t][i];
				}
				model.add(theta[0][i] <= ub);
			}
				



			//Objective
			cout << "Adding objective \n";
			IloExpr obj(env);
			for (int i = 0; i < N; i++) { obj += theta[0][i]; }
			model.add(IloMaximize(env, obj));
			obj.end();

			//Warmstart
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);

			if (!warmstart) {
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						startVal.add((int)data.params["y0"][j]);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (data.params["x0"][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
				}
				for (int i = 0; i < N; i++)
				{
					startVar.add(theta[0][i]);
					startVal.add(0);
				}
			}
			else {
				auto sol = *warmstart;
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						int val = 0;
						if (sol[t][j] > 0) { val = 1; }
						startVal.add(val);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (sol[t][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}

				}
				for (int i = 0; i < N; i++) {
					startVar.add(theta[0][i]);
					IloNum val = 0;
					for (int t = 0; t < T; t++) {
						IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
						for (int r = 0; r < R[i]; r++) {
							if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						}
					}
					startVal.add(val);
				}

			}

			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			break;
		}
			///////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////
		case multicuts::Multi3B1:
		case multicuts::Multi3B2:
		{
			////Upper bound for theta (ensures bounded problem)

			for (int t = 0; t < T; t++) {
				for (int i = 0; i < N; i++) {
					model.add(theta[t][i] <= (IloInt)data.params["Ni"][t][i]);
				}
			}



			//Objective
			cout << "Adding objective \n";
			IloExpr obj(env);
			for (int t = 0; t < T; t++) { 
				for (int i = 0; i < N; i++) {
					obj += theta[t][i];
				}
			}
			model.add(IloMaximize(env, obj));
			obj.end();

			//Warmstart
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);

			if (!warmstart) {
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						startVal.add((int)data.params["y0"][j]);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (data.params["x0"][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
					for (int i = 0; i < N; i++) {
						startVar.add(theta[t][i]);
						startVal.add(0);
					}
				}
			}
			else {
				auto sol = *warmstart;
				for (int t = 0; t < T; t++) {
					for (int j = 0; j < M; j++) {
						startVar.add(y[t][j]);
						int val = 0;
						if (sol[t][j] > 0) { val = 1; }
						startVal.add(val);
						for (int k = 0; k < Mj[j]; k++) {
							int val = 0;
							if (sol[t][j] == k) {
								val = 1;
							}
							startVar.add(x[t][j][k]);
							startVal.add(val);
						}
					}
					for (int i = 0; i < N; i++) {
						startVar.add(theta[t][i]);
						IloNum val = 0;
						IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
						for (int r = 0; r < R[i]; r++) {
							if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						}
						startVal.add(val);
					}
				}

			}

			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			break;
		}
		default:
		{
			break;
		}
		}


		/////////////////////////////////////////////////////////////////////////

		//Link callback
		MulticutCallback cb(data, x, theta, cut_type);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);




		cout << "Model creation time: " << time(NULL) - start << " seconds" << endl;

		//Solve and get results
		IloBool solved = cplex.solve();
		if (solved) {
			cplex.out() << "Solution status:" << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value:" << cplex.getObjValue() << endl;

			ObjectiveValue = cplex.getObjValue();
			SolveTime = cplex.getTime();
			OptimalityGap = cplex.getMIPRelativeGap();


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