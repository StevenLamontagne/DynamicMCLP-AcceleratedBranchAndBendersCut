#include "Model_Multicut.h"
#define EPS 1.0e-6

ILOSTLBEGIN
//typedef IloArray<IloNumArray>    Num2D;
//typedef IloArray<IloNumVarArray>    NumVar2D;
//typedef IloArray<NumVar2D>    NumVar3D;
//typedef IloArray<IloBoolVarArray> BoolVar2D;
//typedef IloArray<BoolVar2D> BoolVar3D;


void Model_Multicut::SetData(const Data& newData)
{
	data = newData; 
	//Define constants for easier reading and writing
	int& T = data.T;
	int& M = data.M;
	int& N = data.N;
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;
	//time_t start;
	//time(&start);

	////Calculate coverage and overlap
	//for (int t = 0; t < T; t++) {
	//	vector<vector<vector<pair<int, int>>>> cover_station1;
	//	for (int j1 = 0; j1 < M; j1++) {
	//		vector<vector<pair<int, int>>> cover_station2;
	//		for (int k1 = 0; k1 < Mj[j1]; k1++) {
	//			vector<pair<int, int>> cover_station3;
	//			cover_station2.push_back(cover_station3);
	//		}
	//		cover_station1.push_back(cover_station2);
	//	}
	//	cover_station.push_back(cover_station1);

	//	for (int i = 0; i < N; i++) {
	//		for (int r = 0; r < R[i]; r++) {
	//			switch (data.P[t][i][r])
	//			{
	//			case triplet::Uncoverable:
	//			case triplet::Precovered: {
	//				continue;
	//				break;
	//			}
	//			default:
	//			{
	//				for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
	//					int j = cover_triplet.first;
	//					int k = cover_triplet.second;
	//					//if (j >= 8|| k >= 4) {
	//					//	cout << "High element found";
	//					//}
	//					cover_station[t][j][k].push_back(make_pair(i, r));
	//				}
	//				break;
	//			}
	//			}

	//		}
	//	}


	//	map<pair<int, int>, map<pair<int, int>, double>> overlap1;
	//	for (int j1 = 0; j1 < M; j1++) {
	//		int total = 0;
	//		vector<pair<int, int>> covered;
	//		for (int k1 = 1; k1 < Mj[j1]; k1++) {
	//			map<pair<int, int>, double> overlap2;
	//			total += cover_station[t][j1][k1].size();
	//			covered.insert(covered.end(), cover_station[t][j1][k1].begin(), cover_station[t][j1][k1].end());
	//			for (int j2 = 0; j2 < M; j2++) {
	//				if (j1 == j2) { continue; }
	//				for (int k2 = 1; k2 < Mj[j2]; k2++) {
	//					int over = 0;
	//					for (pair<int, int> trip : covered) {
	//						int i = trip.first;
	//						int r = trip.second;
	//						if (data.a[t][i][r][j2][k2]) { over += 1; } //not very efficient, since recalculating a ton of stuff. But works for now
	//					}
	//					overlap2[make_pair(j2, k2)] = (double)over / (double)total;
	//				}
	//			}
	//			overlap1[make_pair(j1, k1)] = overlap2;
	//		}
	//	}


	//}

	//cout << "Overlap calculation time: " << time(NULL) - start << " seconds" << endl;
}


void Model_Multicut::Solve(json params)
{
	//Set model parameters
	multicut = params.value("multicut", multicuts::Multi1B1);
	heuristic = params.value("heuristic", useHeuristic::Warmstart);
	nGRASP = params.value("nGRASP", 0);
	verbose = params.value("verbose", false);
	overlap_threshold = params.value("overlap_threshold", 0.1);
	grasp_threshold = params.value("grasp_threshold", 0.85);
	use_trust = params.value("use_trust", false);
	trust_threshold = params.value("trust_threshold", 2.5);

	string outfile = params.value("outfile", "");
	ofstream fout(outfile);

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M = data.M;
	int& N = data.N;
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;
	

	


	time_t start;
	time(&start);
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	if (outfile != "") {
		cplex.setOut(fout);
		cplex.setWarning(fout);
		cplex.setError(fout);
		cplex.out() << "Log sent to file" << endl;
	}
	else {
		cout << "No log file provided" << endl;
	}
	try {
		//Create model and set parameters
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200); 
		cplex.setParam(IloCplex::Param::Preprocessing::Linear, 0);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);
		if (!verbose) { 
			cplex.setParam(IloCplex::Param::MIP::Display, 0); 
			cplex.setParam(IloCplex::Param::Tune::Display, 0);
			cplex.setParam(IloCplex::Param::Simplex::Display, 0);
			cplex.setParam(IloCplex::Param::Sifting::Display, 0);
			cplex.setParam(IloCplex::Param::ParamDisplay, 0);
			cplex.setParam(IloCplex::Param::Network::Display, 0);
			cplex.setParam(IloCplex::Param::Conflict::Display, 0);
		}
		else { 
			//cplex.setParam(IloCplex::Param::MIP::Display, 1);
			cplex.setParam(IloCplex::Param::MIP::Interval, 10000);
			cplex.setParam(IloCplex::Param::Tune::Display, 0);
			cplex.setParam(IloCplex::Param::Simplex::Display, 0);
			cplex.setParam(IloCplex::Param::Sifting::Display, 0);
			cplex.setParam(IloCplex::Param::ParamDisplay, 0);
			cplex.setParam(IloCplex::Param::Network::Display, 0);
			cplex.setParam(IloCplex::Param::Conflict::Display, 0);
		}

		//Create variables
		BoolVar3D x(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = BoolVar2D(env, M);
			for (int j = 0; j < M; j++) {
				x[t][j] = IloBoolVarArray(env, Mj[j]);
			}
		}


		NumVar2D theta(env, T);
		for (int t = 0; t < T; t++) {
			theta[t] = IloNumVarArray(env, N, 0.0, IloInfinity); //NumVar initialisation requires the bounds, or accessing it will crash!
		}

		IloNumVar theta_obj(env);


		//Constraints
		////Budget, year 0
		IloExpr budget0(env);
		for (int j = 0; j < M; j++) {
			for (int k = 1; k < Mj[j]; k++) {
				budget0 += (double)data.params["c"][0][j][k] * (x[0][j][k] - (int)data.params["x0"][j][k]);
			}
		}
		model.add(budget0 <= (double)data.params["B"][0]);
		budget0.end();

		////Budget, year 1+
		for (int t = 1; t < T; t++) {
			IloExpr budget(env);
			for (int j = 0; j < M; j++) {
				for (int k = 1; k < Mj[j]; k++) {
					budget += (double)data.params["c"][t][j][k] * (x[t][j][k] - x[t-1][j][k]);
				}
			}
			model.add(budget <= (double)data.params["B"][t]);
			budget.end();
		}

		////Can't remove outlets, year 0
		for (int j = 0; j < M; j++) {
			for (int k = 1; k < Mj[j]; k++) {
				model.add(x[0][j][k] >= (int)data.params["x0"][j][k]);
			}
		}

		////Can't remove outlets, year 1+
		for (int t = 1; t < T; t++) {
			for (int j = 0; j < M; j++) {
				for (int k = 1; k < Mj[j]; k++) {
					model.add(x[t][j][k] >= x[t-1][j][k]);
				}
			}
		}

		////At least k outlets, year 0+
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				for (int k = 1; k < Mj[j]; k++) {
					model.add(x[t][j][k] <= x[t][j][k - 1]);
				}
			}
		}

		///Redundant, but prevents crash from unextracted variables
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				for (int k = 0; k < Mj[j]; k++) {
					model.add(x[t][j][k] >= 0);
				}		
			}
		}


		//Set objective via proxy. Cut methods define theta_obj via theta variables
		model.add(IloMaximize(env, theta_obj));



		
		//Create warmstart solution (defaults to zero solution if not greedy)
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		vector<vector<int>> sol;

		if (heuristic == useHeuristic::Warmstart || heuristic == useHeuristic::WarmstartAndPostSimple || heuristic == useHeuristic::WarmstartAndPostGreedy) {
			{Greedy G;
			G.SetData(data);
			G.Solve();
			sol = G.Solution;
			ObjectiveValue = G.SolutionQuality;
			//for (int t = 0; t < T; t++) {
			//	vector<int> sol1;
			//	for (int j = 0; j < M; j++) {
			//		if (G.Solution[t][j] >= j) { sol1.push_back(1); }
			//		else { sol1.push_back(0); }
			//	}
			//	sol.push_back(sol1);
			//}
			}
		}
		else {
			for (int t = 0; t < T; t++) {
				vector<int> sol1;
				for (int j = 0; j < M; j++) {
					sol1.push_back(0);
				}
				sol.push_back(sol1);
			}
		}

		//Cut specific constraints and objective
		/////////////////////////////////////////////////////////////////////////
		switch (multicut)
		{
			///////////////////////////////////////////////////////////
		case multicuts::SingleB1:
		case multicuts::SingleB2:
		case multicuts::SinglePO1:
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
			if (verbose) { cplex.out() << "Adding objective \n"; }
			IloExpr obj(env);
			obj += theta[0][0];
			model.add(theta_obj == obj);
			obj.end();

			for (int t = 0; t < T; t++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < Mj[j]; k++) {
						int val = 0;
						if (sol[t][j] >= k) {
							val = 1;
						}
						startVar.add(x[t][j][k]);
						startVal.add(val);
					}
				}
			}
			startVar.add(theta[0][0]);
			startVal.add(data.SolutionQuality(sol));

			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			break;
		}
			///////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////
		case multicuts::Multi1SB1:
		case multicuts::Multi1B1:
		case multicuts::Multi1B2:
		case multicuts::Multi1PO1:
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
			if (verbose) { cplex.out() << "Adding objective \n"; }
			IloExpr obj(env);
			for (int t = 0; t < T; t++) { obj += theta[t][0]; }
			model.add(theta_obj == obj);
			obj.end();

			//Warmstart
			for (int t = 0; t < T; t++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < Mj[j]; k++) {
						int val = 0;
						if (sol[t][j] >= k) {
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
						if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
						if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						else {
							bool covered = 0;
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]){
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
			if (verbose) { cplex.out() << "Adding objective \n"; }
			IloExpr obj(env);
			for (int i = 0; i < N; i++) { obj += theta[0][i]; }
			model.add(theta_obj == obj);
			obj.end();

			//Warmstart
			for (int t = 0; t < T; t++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < Mj[j]; k++) {
						int val = 0;
						if (sol[t][j] >= k) {
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
						if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
						if (data.P[t][i][r] == triplet::Precovered) { val += weight; }
						else {
							bool covered = 0;
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]){
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
		case multicuts::Multi3SB2:
		case multicuts::Multi3B1:
		case multicuts::Multi3B2:
		case multicuts::Multi3PO1:
		{
			////Upper bound for theta (ensures bounded problem)

			for (int t = 0; t < T; t++) {
				for (int i = 0; i < N; i++) {
					model.add(theta[t][i] <= (IloInt)data.params["Ni"][t][i]);
				}
			}



			//Objective
			if (verbose) { cplex.out() << "Adding objective \n"; }
			IloExpr obj(env);
			for (int t = 0; t < T; t++) { 
				for (int i = 0; i < N; i++) {
					obj += theta[t][i];
				}
			}
			model.add(theta_obj == obj);
			obj.end();

			//Warmstart
			for (int t = 0; t < T; t++) {
				for (int j = 0; j < M; j++) {
					for (int k = 0; k < Mj[j]; k++) {
						int val = 0;
						if (sol[t][j] >= k) {
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]){
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
			cout << "You forgot to add the warmstart section for this option, dummy" << endl;
			throw exception();
			break;
		}
		}

		//Add GRASP cuts
		if (nGRASP > 0) {
			if (verbose) { cplex.out() << "Adding GRASP cuts" << endl; }
			int success = 0;
			int failure = 0;
			for (int temp = 0; temp < nGRASP; temp++) {
				bool s = GRASPCut(env, model, x, theta);
				if (s) { success += 1; }
				else {
					failure += 1;
					if (failure >= success + 10) { break; }
				}
			}

		}

		/////////////////////////////////////////////////////////////////////////

		//Link callback
		MulticutCallback cb(data, x, theta, multicut, heuristic);
		cb.threshold = overlap_threshold;
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);


		bool skip_solve = false;
		//Trust region
		if (use_trust) {
			IloConstraintArray trust(env);

			//Note: Trust distance definitions depend only on the warmstart solution, so never need to be updated
			IloExpr trust_expr(env);
			NumVar3D delta(env, T);
			//Num3D I_tilde(env, T);
			for (int t = 0; t < T; t++) {
				delta[t] = NumVar2D(env, M);
				//I_tilde[t] = Num2D(env, N);
				for (int j = 0; j < M; j++) {
					delta[t][j] = IloNumVarArray(env, Mj[j], 0.0, IloInfinity); //NumVar initialisation requires the bounds, or accessing it will crash!
					for (int k = 1; k < Mj[j]; k++) {
						IloExpr expr(env);
						expr -= delta[t][j][k];
						for (int t_bar = 0; t_bar < t; t_bar++) { expr += delta[t_bar][j][k]; } //We add the previous delta variables to account for outlet preservation constraints
						if (sol[t][j] >= k) { expr += 1 - x[t][j][k];	}
						else { expr += x[t][j][k] ;}
						model.add(expr == 0);
						expr.end();
						trust_expr += delta[t][j][k];
					}
				}

				//for (int i = 0; i < N; i++) {
				//	I_tilde[t][i] = IloNumArray(env, R[i]);
				//	for (int r = 0; r < R[i]; r++) {
				//		IloInt val = 0;
				//		for (auto cover_triplet : data.cover_triplet[t][i][r]) {
				//			int j = cover_triplet.first;
				//			int k0 = cover_triplet.second;
				//			if (sol[t][j] >= k0) { val += 1; }
				//		}
				//		I_tilde[t][i][r] = val;
				//	}
				//}
			}


			//for (int t = 0; t < T; t++) {
			//	int temp = 0;
			//	for (int j = 0; j < M; j++) {
			//		if (sol[t][j] < 1) { continue; }
			//		Num3D I_hat(env);
			//		I_hat = I_tilde;
			//		IloNum obj = ObjectiveValue;

			//		for (int t2 = t; t2 < T; t2++) {
			//			int current = sol[t2][j];
			//			for (auto cover_station : data.cover_station[t2][j][current]) {
			//				int i = cover_station.first;
			//				int r = cover_station.second;
			//				I_hat[t2][i][r] -= 1;
			//			}
			//			int temp = 0;
			//		}
			//		I_hat.end();
			//	}
			//}

			{
				IloConstraint cons = trust_expr <= trust_threshold;
				model.add(cons);
				trust.add(cons);
				//cons.end();
			}

			IloBool trust_solved = cplex.solve();
			if (trust_solved) {
				if (verbose) {
					cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
					cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
					cplex.out() << "Number of nodes: " << cplex.getNnodes() << endl;
					cplex.out() << "Trust region solve time: " << cplex.getTime() << endl;
				}
				ObjectiveValue = cplex.getObjValue();
				SolveTime = cplex.getTime();
				TotalTime = time(NULL) - start;
				OptimalityGap = cplex.getMIPRelativeGap();
				GetSolution(cplex, x);
				stats = cb.stats;
				stats["nNodes"] = cplex.getNnodes();
				stats["TrustSolveTime(DivideBy100)"] = (int) (100*cplex.getTime());
			}
			//I_tilde.end();
			IloNum timelimit = 7200 - cplex.getTime();
			if (timelimit >EPS){ cplex.setParam(IloCplex::Param::TimeLimit, timelimit); }
			else { skip_solve = true; }
			//for (int i = 0; i < trust.getSize(); i ++) {
			//	model.remove(trust[i]);
			//}
			cplex.clearUserCuts();
			cplex.deleteMIPStarts(0, cplex.getNMIPStarts());
			{
				model.remove(trust);
				IloConstraint cons = trust_expr >= trust_threshold + 1;
				model.add(cons);
				trust.add(cons);
			}

			IloConstraint cons;
			IloExpr obj(env);
			switch (multicut)
			{

			case multicuts::SingleB1:
			case multicuts::SingleB2:
			case multicuts::SinglePO1:
				cons = theta[0][0] >= ObjectiveValue;
				model.add(cons);
				break;
			case multicuts::Multi1SB1:
			case multicuts::Multi1B1:
			case multicuts::Multi1B2:
			case multicuts::Multi1PO1:
				for (int t = 0; t < T; t++) { obj += theta[t][0]; }
				cons = obj >= ObjectiveValue;
				model.add(cons);
				break;
			case multicuts::Multi2B1:
			case multicuts::Multi2B2:
				for (int i = 0; i < N; i++) { obj += theta[0][i]; }
				cons = obj >= ObjectiveValue;
				model.add(cons);
				break;
			case multicuts::Multi3SB2:
			case multicuts::Multi3B1:
			case multicuts::Multi3B2:
			case multicuts::Multi3PO1:
				for (int t = 0; t < T; t++) { for (int i = 0; i < N; i++) { obj += theta[t][i]; } }
				cons = obj >= ObjectiveValue;
				model.add(cons);
				break;
			default:
				break;
			}
		}

		if (verbose) { cplex.out() << "Model creation time: " << time(NULL) - start << " seconds" << endl; }

		//Solve and get results
		if (!skip_solve) {
			IloBool solved = cplex.solve();
			if (solved) {
				if (verbose) {
					cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
					cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
					cplex.out() << "Number of nodes: " << cplex.getNnodes() << endl;
					cplex.out() << "\n" << endl;
				}

				if (use_trust) {
					for (pair<string, int> res : cb.stats) {
						string category = res.first;
						int value = res.second;
						stats[category] = stats[category] + value;
					}
					stats["nNodes"] += cplex.getNnodes();
				}
				else {
					stats = cb.stats;
					stats["nNodes"] = cplex.getNnodes();
				}

				ObjectiveValue = cplex.getObjValue();
				SolveTime = cplex.getTime();
				if (use_trust) { SolveTime += ((double)stats["TrustSolveTime(DivideBy100)"]) / 100; }
				TotalTime = time(NULL) - start;
				OptimalityGap = max((double)OptimalityGap, cplex.getMIPRelativeGap());
				GetSolution(cplex, x);

			}
			else if (!use_trust) {
				ObjectiveValue = -1;
				SolveTime = -1;
				OptimalityGap = -1;
				stats = cb.stats;
				stats["nNodes"] = cplex.getNnodes();
			}
		}
	}
	catch (IloException & e) {
		cplex.out() << "Exception: " << e << endl;
		ObjectiveValue = -1;
		SolveTime = -1;
		OptimalityGap = -1;
	}
	cplex.clearUserCuts();
	cplex.end();
	model.end();
	env.end();
	fout.close();

}

void Model_Multicut::GetSolution(IloCplex& cplex, BoolVar3D& x)
{
	vector<vector<int>> temp;
	for (int t = 0; t < data.T; t++) {
		vector<int> Solution1;
		for (int j = 0; j < data.M; j++) {
			int total = 0;
			for (int k = 1; k < data.Mj[j]; k++) {
				if (cplex.getValue(x[t][j][k]) > 1 - EPS) {
					total = k;
				}
			}
			Solution1.push_back(total);
		}
		temp.push_back(Solution1);
	}
	Solution = temp;
}



/*
Determine greedy randomised solution and add associated cut
*/
bool Model_Multicut::GRASPCut(IloEnv& env, IloModel& model, BoolVar3D& x, IloArray<IloNumVarArray>& theta)
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
					if (I_tilde[t][i][r] >=1) { continue; }
					for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
						int j = cover_triplet.first;
						int k0 = cover_triplet.second;
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
				double threshold = grasp_threshold * maxVal;
				vector<int> candidates;
				for (int j = 0; j < data.M; j++) {
					if (StationTotals[j] >= threshold) { candidates.push_back(j);  }
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
							if ( (I_tilde[tBar][i][r] <=1) && (data.a[tBar][i][r][j_star][Sol[tBar][j_star]])) {
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
							int j = data.cover_triplet[t][i][r][0].first;
							int k0 = data.cover_triplet[t][i][r][0].second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							int j = data.cover_triplet[t][i][r][0].first;
							int k0 = data.cover_triplet[t][i][r][0].second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							int j = data.cover_triplet[t][i][r][0].first;
							int k0 = data.cover_triplet[t][i][r][0].second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							int j = data.cover_triplet[t][i][r][0].first;
							int k0 = data.cover_triplet[t][i][r][0].second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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
							for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
								int j = cover_triplet.first;
								int k0 = cover_triplet.second;
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






