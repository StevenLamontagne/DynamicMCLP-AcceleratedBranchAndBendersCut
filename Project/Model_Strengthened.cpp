#include "Model_Strengthened.h"
#define EPS 1.0e-6

ILOSTLBEGIN
typedef IloArray<IloNumArray>    Float2D;
typedef IloArray<IloNumVarArray>    FloatVar2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;


void Model_Strengthened::Solve(multicuts cut_type, useHeuristic heuristic, int nGRASP)
{
	multicut = cut_type;
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
		cplex.setParam(IloCplex::Param::Preprocessing::Linear , 0);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);

		//Create variables
		BoolVar3D x(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = BoolVar2D(env, M);
			for (int j = 0; j < M; j++) {
				x[t][j] = IloBoolVarArray(env, Mj[j]);
			}
		}

		BoolVar3D z(env, T);
		for (int t = 0; t < T; t++) {
			z[t] = BoolVar2D(env, M);
			for (int j = 0; j < M; j++) {
				z[t][j] = IloBoolVarArray(env, Mj[j]);
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
		cout << "Adding constraints" << endl;

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
				open.end();
			}
		}


		//// Set z
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				////Set based on higher values
				//for (int k = 1; k < Mj[j]; k++) {
				//	IloExpr temp(env);
				//	for (int k1 = k; k1 < Mj[j]; k1++) {
				//		temp += x[t][j][k1];
				//	}
				//	model.add(z[t][j][k] == temp);
				//	temp.end();

				//	model.add(z[t][j][k] <= z[t][j][1]); //Redudant in theory, but prevents the solver from trying shenanigans for the budget
				//}

				////Set based on number of outlets
				//IloExpr temp(env);
				//for (int k = 1; k < Mj[j]; k++) {
				//	temp += k * x[t][j][k];
				//}

				//for (int k = 1; k < Mj[j]; k++) {
				//	model.add(k*z[t][j][k] <= temp);
				//	model.add(z[t][j][k] <= z[t][j][1]); //Redudant in theory, but prevents the solver from trying shenanigans for the budget
				//}
				//temp.end();


				////Set based on number of outlets (better)
				for (int k = 1; k < Mj[j]; k++) {
					IloExpr temp(env);
					for (int kBar = k; kBar < Mj[j]; kBar++) {
						temp += kBar * x[t][j][kBar];
					}
					model.add(k*z[t][j][k] <= temp);
					model.add(z[t][j][k] <= z[t][j][1]); //Redudant in theory, but prevents the solver from trying shenanigans for the budget
					temp.end();
				}
			}
		}

		////Pay one-time cost
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				model.add(z[t][j][1] <= y[t][j]);
			}
		}

		///Redundant, but prevents crash from unextracted variables
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				model.add(theta[t][i] >= 0);
			}
			for (int j = 0; j < M; j++) {
				model.add(y[t][j] >= 0);
				for (int k = 0; k < Mj[j]; k++) {
					model.add(x[t][j][k] >= 0);
					model.add(z[t][j][k] >= 0);
				}
			}
		}

		////Sets x variables as SOS1
		//for (int t= 0; t < T; t++) {
		//	for (int j = 0; j < M; j++) {
		//		IloNumVarArray sosvar(env, Mj[j]);
		//		IloNumArray    sosval(env, Mj[j]);
		//		for (int k = 0; k < Mj[j]; k++) {
		//			sosvar[k] = x[t][j][k];
		//			sosval[k] = k;
		//		}
		//		model.add(IloSOS1(env, sosvar, sosval));
		//	}
		//}





		//Cut specific constraints and objective
		////Create warmstart solution
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		vector<vector<int>> sol;
		for (int t = 0; t < T; t++) {
				vector<int> sol1;
				for (int j = 0; j < M; j++) {
					sol1.push_back(0);
				}
				sol.push_back(sol1);
		}

		/////////////////////////////////////////////////////////////////////////
		switch (cut_type)
		{
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
			for (int t = 0; t < T; t++) { 
				obj += theta[t][0]; 
			}
			model.add(IloMaximize(env, obj));
			obj.end();

			//Warmstart
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

					IloNum outlets = 0.0;
					for (int k = 0; k < Mj[j]; k++) {
						startVar.add(z[t][j][k]);
						if (sol[t][j] >= k) { startVal.add(1); }
						else { startVal.add(0); }

					}
				}



				startVar.add(theta[t][0]);
				IloNum val = 0;
				for (int i = 0; i < N; i++) {
					IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
					for (int r = 0; r < R[i]; r++) {
						if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
						if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						else {
							bool covered = 0;
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								if (sol[t][j] >= k0) { covered = 1; }
							}
							if (covered) { val += weight; }
						}
					}
				}
				startVal.add(val);
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
						if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
						if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						else {
							bool covered = 0;
							vector<pair<int, int>> cover = data.cover[t][i][r];
							for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
								int j = cover[jBar].first;
								int k0 = cover[jBar].second;
								if (sol[t][j] >= k0) { covered = 1; }
							}
							if (covered) { val += weight; }
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
		default:
		{
			break;
		}
		}

		/////////////////////////////////////////////////////////////////////////

		//Link callback
		MulticutCallback_Strengthened cb(data, z, y, theta, cut_type, heuristic);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);




		cout << "Model creation time: " << time(NULL) - start << " seconds" << endl;

		//Solve and get results
		IloBool solved = cplex.solve();
		if (solved) {
			cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
			//cplex.out() << "Repairs attempted: " << cb.totalCounter << endl;
			//cplex.out() << "Repairs succeeded: " << cb.repairCounter << endl;


			ObjectiveValue = cplex.getObjValue();
			SolveTime = cplex.getTime();
			TotalTime = time(NULL) - start;
			OptimalityGap = cplex.getMIPRelativeGap();
			if (cb.totalCounter == 0) { RepairRatio = -1; }
			else { RepairRatio = (double)cb.repairCounter / (double)cb.totalCounter; }


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

bool Model_Strengthened::GreedyFill(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage)
{
	////Initialise solution and budget
	//for (int t = 0; t < T; t++) {
	//	for (int j = 0; j < M; j++) {
	//		for (int k = 1; k < Mj[j]; k++) {
	//			if (oldSolution[t][j][k] > EPS) {
	//				newSolution[t][j] = k - 1 + oldSolution[t][j][k];
	//				break;
	//			}
	//		}

	//		double previous = 0.0;
	//		if (t == 0) { previous = data.params["x0"][j]; }
	//		else { previous = newSolution[t-1][j]; }

	//		double diff = newSolution[t][j] - previous;
	//		if ((diff > EPS) && (previous < 1- EPS)) {
	//			Budget[t] -= diff * data.params["f"][t][j];
	//		}
	//		if (diff > EPS) {
	//			Budget[t] -= diff * data.params["c"][t][j];
	//		}
	//	}
	//}
	//if (sum(Budget) < EPS) { return false; }

	//Check if solution is fractional and, if so, attempt to repair it


	//Exit unsuccessfully if solution remains fractional


	//Otherwise, begin greedy fill

	//vector<vector<vector<double>>> coverage;
	//double SolutionQuality = 0.0;

	////Initialisation
	//for (int t = 0; t < T; t++) {
	//	//Initialise coverage of triplets
	//	vector<vector<double>> cover1;
	//	for (int i = 0; i < N; i++) {
	//		vector<double> cover2;
	//		double weight = (double)data.params["Ni"][t][i] / (double)R[i];
	//		for (int r = 0; r < R[i]; r++) {
	//			bool val = 0;
	//			switch (data.P[t][i][r])
	//			{
	//			case triplet::Uncoverable:
	//				val = 1; //Set to 1 to skip trying to cover later
	//				break;
	//			case triplet::Precovered:
	//				val = 1;
	//				SolutionQuality += weight;
	//				break;
	//			default:
	//				break;
	//			}
	//			cover2.push_back(val);
	//		}
	//		cover1.push_back(cover2);
	//	}
	//	coverage.push_back(cover1);
	//}

	//Initialise vectors for holding information for each station
	vector<double> StationTotals;
	vector<double> Costs;

	for (int j = 0; j < data.M; j++) {
		StationTotals.push_back(0.0);
		Costs.push_back(0.0);
	}

	bool foundImprovement;
	//Greedy optimisation
	for (int t = 0; t < data.T; t++) {
		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			for (int j = 0; j < data.M; j++) {
				StationTotals[j] = 0.0;
				Costs[j] = 0.0;
			}

			//Check for coverage by adding outlet
			for (int i = 0; i < data.N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (coverage[t][i][r]) { continue; }
					vector<pair<int, int>> cover = data.cover[t][i][r];
					for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
						int j = cover[jBar].first;
						int k0 = cover[jBar].second;
						if (Sol[t][j] + 1 >= k0) {
							StationTotals[j] += weight;
						}
					}
				}
			}

			//Reset coverage of stations to 0 if can't add new outlet
			for (int j = 0; j < data.M; j++) {
				//Maximum number of outlets already placed
				if (Sol[t][j] >= data.Mj[j] - 1) {
					StationTotals[j] = 0.0;
					continue;
				}

				//Budget insufficient for new outlet
				Costs[j] += (double)data.params["c"][t][j];
				if (Sol[t][j] == 0) { Costs[j] += (double)data.params["f"][t][j]; }
				if (Costs[j] > Budget[t]) {
					StationTotals[j] = 0.0;
				}
				//cout << "Coverage of station " << j << ": " << StationTotals[j] << endl;
			}


			int j_star = argmax(StationTotals);
			double z_star = StationTotals[j_star];
			if (z_star > 0) {
				//Update budget
				Budget[t] -= Costs[j_star];

				for (int tBar = t; tBar < data.T; tBar++) {
					//Update solution for current (and all future) years
					Sol[tBar][j_star] += 1;

					//Update coverage of triplets
					for (int i = 0; i < data.N; i++) {
						//double weight = (double)data.params["Ni"][t][i] / (double)R[i];
						for (int r = 0; r < data.R[i]; r++) {
							if ((!coverage[tBar][i][r]) && (data.aBar[tBar][i][r][j_star][Sol[tBar][j_star]])) {
								coverage[tBar][i][r] = 1;

							}
						}
					}
				}
				foundImprovement = 1;
			}


		} while (foundImprovement);
	}

	return false;
}



vector<int> Model_Strengthened::GetFractional(vector<vector<double>> Sol)
{
	vector<int> fracList;
	for (int j = 0; j < data.M; j++) {
		bool isFrac = false;
		for (int t = 0; t < data.T; t++) {
			double current;
			double frac = modf(Sol[t][j], &current);
			if (frac > EPS) { isFrac = true; }
		}
		if (isFrac) { fracList.push_back(j); }
	}
	return fracList;
}

bool Model_Strengthened::GRASPCut(IloEnv& env, IloModel& model, BoolVar3D& x, IloArray<IloNumVarArray>& theta)
{
	//Initialising baseline information
	vector<double> Budget = data.params["B"];

	IloArray<IloArray<IloIntArray>> I_tilde(env, data.T);
	for (int t = 0; t < data.T; t++) {
		I_tilde[t] = IloArray<IloIntArray>(env, data.N);
		for (int i = 0; i < data.N; i++) {
			I_tilde[t][i] = IloIntArray(env, data.R[i]);
		}
	}

	IloArray<IloIntArray> Sol(env, data.T);
	for (int t = 0; t < data.T; t++) {
		Sol[t] = IloIntArray(env, data.M);
		for (int j = 0; j < data.M; j++) {
			Sol[t][j] = data.params["x0"][j];
		}
	}


	//Initialise vectors for holding information for each station
	vector<double> StationTotals;
	vector<double> Costs;
	for (int j = 0; j < data.M; j++) {
		StationTotals.push_back(0.0);
		Costs.push_back(0.0);
	}

	//GRASP optimisation
	for (int t = 0; t < data.T; t++) {
		bool foundImprovement;

		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			for (int j = 0; j < data.M; j++) {
				StationTotals[j] = 0.0;
				Costs[j] = 0.0;
			}

			//Check for coverage by adding outlet
			for (int i = 0; i < data.N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (I_tilde[t][i][r] >= 1) { continue; }
					vector<pair<int, int>> cover = data.cover[t][i][r];
					for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
						int j = cover[jBar].first;
						int k0 = cover[jBar].second;
						if (Sol[t][j] + 1 >= k0) {
							StationTotals[j] += weight;
						}
					}
				}
			}

			//Reset coverage of stations to 0 if can't add new outlet
			for (int j = 0; j < data.M; j++) {
				//Maximum number of outlets already placed
				if (Sol[t][j] >= data.Mj[j] - 1) {
					StationTotals[j] = 0.0;
					continue;
				}

				//Budget insufficient for new outlet
				Costs[j] += (double)data.params["c"][t][j];
				if (Sol[t][j] == 0) { Costs[j] += (double)data.params["f"][t][j]; }
				if (Costs[j] > Budget[t]) {
					StationTotals[j] = 0.0;
				}
			}

			double maxVal = *max_element(StationTotals.begin(), StationTotals.end());
			if (maxVal > 0) {
				double threshold = 0.85 * maxVal;
				vector<int> candidates;
				for (int j = 0; j < data.M; j++) {
					if (StationTotals[j] >= threshold) { candidates.push_back(j); }
				}
				int randIndex = rand() % candidates.size();
				int j_star = candidates[randIndex];

				//Update budget
				Budget[t] -= Costs[j_star];

				for (int tBar = t; tBar < data.T; tBar++) {
					//Update solution for current (and all future) years
					Sol[tBar][j_star] += 1;

					//Update coverage of triplets
					for (int i = 0; i < data.N; i++) {
						for (int r = 0; r < data.R[i]; r++) {
							//We cap the I_tilde at 2, since all values > 1 will result in the same contribution for cuts
							if ((I_tilde[tBar][i][r] <= 1) && (data.aBar[tBar][i][r][j_star][Sol[tBar][j_star]])) {
								I_tilde[tBar][i][r] += 1;
							}
						}
					}
				}
				foundImprovement = 1;
			}


		} while (foundImprovement);
	}


	//Create cuts

	switch (multicut)
	{
		//////////////////////////////////////////////////////////////////////////////////////////
		//Single B1 cut
	case multicuts::SingleB1:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

		for (int t = 0; t < data.T; t++) {
			for (int i = 0; i < data.N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < data.Mj[j]; k++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		model.add(lhs >= 0).end();
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B2 cut
	case multicuts::SingleB2:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

		for (int t = 0; t < data.T; t++) {
			for (int i = 0; i < data.N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		model.add(lhs >= 0).end();
	}
	//////////////////////////////////////////////////////////////////////////////////////////
		//Multi-cut (by year), B1
	case multicuts::Multi1B1: {
		for (int t = 0; t < data.T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < data.N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < data.Mj[j]; k++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			model.add(lhs >= 0).end();
		}

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year), B2
	case multicuts::Multi1B2: {
		for (int t = 0; t < data.T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < data.N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			model.add(lhs >= 0).end();
		}

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by user class), B1
	case multicuts::Multi2B1: {
		for (int i = 0; i < data.N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < data.T; t++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < data.Mj[j]; k++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			model.add(lhs >= 0).end();
		}

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by user class), B2
	case multicuts::Multi2B2: {
		for (int i = 0; i < data.N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < data.T; t++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			model.add(lhs >= 0).end();
		}

		break; //End of switch statement
	}
							//////////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year and user class), B1
	case multicuts::Multi3B1: {
		for (int t = 0; t < data.T; t++) {
			for (int i = 0; i < data.N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < data.Mj[j]; k++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				model.add(lhs >= 0).end();
			}
		}

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year and user class), B2
	case multicuts::Multi3B2: {
		for (int t = 0; t < data.T; t++) {
			for (int i = 0; i < data.N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
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
								for (int k = k0; k < data.Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				model.add(lhs >= 0).end();
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
	return true;
}


bool Model_Strengthened::GreedyRepair(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage)
{
	vector<int> frac = GetFractional(Sol);
	if (frac.size() == 0) {
		//Can skip the greedy fill phase: no improvements are possible
		if (sum(Budget) < EPS) { return false; }
		//Repair unnecessary, but greedy fill phase may find improvements
		else { return true; }
	}


	vector<vector<double>> newSolution = Sol;

	//Refund budget for fractional solutions and set values to 0
	for (int j : frac) {
		for (int t = 0; t < data.T; t++) {
			double previous = 0.0;
			if (t == 0) { previous = data.params["x0"][j]; }
			else { previous = Sol[t - 1][j]; }

			double delta = Sol[t][j] - previous;
			if (previous < 1) { Budget[t] += delta * (double)data.params["f"][t][j]; }
			Budget[t] += delta * (double)data.params["c"][t][j];

			newSolution[t][j] = 0;
		}

	}



	for (int t = 0; t < data.T; t++) {
		//Set all stations to lower bounds, adjust budget
		for (int j : frac) {
			double previous = 0.0;
			if (t == 0) { previous = data.params["x0"][j]; }
			else { previous = newSolution[t - 1][j]; }

			int lb = max(previous, floor(Sol[t][j]));

			if (lb == previous) { continue; }
			int delta = lb - previous;
			if (delta > 0) {
				if (previous == 0) { Budget[t] -= (double)data.params["f"][t][j]; }
				Budget[t] -= (double)data.params["c"][t][j] * delta;
			}
		}

		//If budget is negative, immediately exit and skip greedy fill
		if (Budget[t] < EPS) { return false; }


		vector<double> StationTotals;
		vector<double> Costs;
		for (long unsigned int jBar = 0; jBar < frac.size(); jBar++) {
			StationTotals.push_back(0.0);
			Costs.push_back(0.0);
		}

		bool foundImprovement;
		do
		{
			foundImprovement = false;
			for (long unsigned int jBar = 0; jBar < frac.size(); jBar++) {
				StationTotals[jBar] = 0.0;
				Costs[jBar] = 0.0;
			}


			//Check for coverage from new stations
			for (int i = 0; i < data.N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (coverage[t][i][r]) { continue; }
					for (long unsigned int jBar = 0; jBar < frac.size(); jBar++) {
						int j = frac[jBar];
						if (data.aBar[t][i][r][j][newSolution[t][j] + 1]) {
							StationTotals[jBar] += weight;
						}
					}
				}
			}

			//Reset coverage of stations to 0 if can't add new outlet
			for (long unsigned int jBar = 0; jBar < frac.size(); jBar++) {
				int j = frac[jBar];
				//Ceiling of current value reached
				if (newSolution[t][j] >= ceil(Sol[t][j])) {
					StationTotals[jBar] = 0.0;
					continue;
				}

				//Budget insufficient for new outlet
				Costs[jBar] += (double)data.params["c"][t][j];
				if (newSolution[t][j] == 0) { Costs[jBar] += (double)data.params["f"][t][j]; }
				if (Costs[jBar] > Budget[t]) {
					StationTotals[jBar] = 0.0;
				}
			}

			int jBar_star = argmax(StationTotals);
			double zBar_star = StationTotals[jBar_star];
			if (zBar_star > 0) {
				int j_star = frac[jBar_star];

				//Update budget
				Budget[t] -= Costs[jBar_star];

				for (int tBar = t; tBar < data.T; tBar++) {
					//Update solution for current (and all future) years
					newSolution[tBar][j_star] += 1;

					//Update coverage of triplets
					for (int i = 0; i < data.N; i++) {
						for (int r = 0; r < data.R[i]; r++) {
							if ((!coverage[tBar][i][r]) && (data.aBar[tBar][i][r][j_star][newSolution[tBar][j_star]])) {
								coverage[tBar][i][r] = 1;
							}
						}
					}
				}
				foundImprovement = 1;

			}



		} while (foundImprovement);
	}

	Sol = newSolution;
	return true;
}



