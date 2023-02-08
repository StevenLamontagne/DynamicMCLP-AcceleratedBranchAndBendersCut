#include "Model_Improved.h"

#define EPS 1.0e-6

ILOSTLBEGIN
//typedef IloArray<IloNumArray>    Num2D;
//typedef IloArray<IloNumVarArray>    NumVar2D;
//typedef IloArray<NumVar2D>    NumVar3D;
//typedef IloArray<IloBoolVarArray> BoolVar2D;
//typedef IloArray<BoolVar2D> BoolVar3D;


void Model_Improved::SetData(const Data_Improved& newData)
{
	data = newData;

	////Define constants for easier reading and writing
	//int& T = data.T;
	//int& M_bar = data.M_bar;
	//vector<int>& P = data.P;

	//int M = data.params["M"];
	//vector<int> Mj = data.params["Mj"];
	//map<pair<int, int>, int> coordConverter;
	//for (int j_bar = 0; j_bar < M_bar; j_bar++) { coordConverter[make_pair(data.params["station_coord"][j_bar][0], data.params["station_coord"][j_bar][1])] = j_bar; }

}


void Model_Improved::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);
	overlap_threshold = params.value("overlap_threshold", 0.1);
	use_trust = params.value("use_trust", false);
	trust_threshold = params.value("trust_threshold", 2.0);

	string outfile = params.value("outfile", "");
	ofstream fout(outfile);

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	//int M = data.params["M"];
	//vector<int> Mj;
	//for (int val : data.params["Mj"]) { Mj.push_back(val); }





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
		BoolVar2D x(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = IloBoolVarArray(env, M_bar);
		}


		IloNumVarArray  theta(env, T, 0.0, IloInfinity);

		IloNumVar theta_obj(env);


		//Constraints
		////Budget, year 0
		IloExpr budget0(env);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			budget0 += (double)data.params["c"][0][j][k] * (x[0][j_bar] - (int)data.params["x0"][j][k]);
		}
		model.add(budget0 <= (double)data.params["B"][0]);
		budget0.end();

		////Budget, year 1+
		for (int t = 1; t < T; t++) {
			IloExpr budget(env);
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				budget += (double)data.params["c"][t][j][k] * (x[t][j_bar] - x[t - 1][j_bar]);
			}
			model.add(budget <= (double)data.params["B"][t]);
			budget.end();
		}

		////Can't remove outlets, year 0
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			model.add(x[0][j_bar] >= (int)data.params["x0"][j][k]);
		}

		////Can't remove outlets, year 1+
		for (int t = 1; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				model.add(x[t][j_bar] >= x[t - 1][j_bar]);
			}
		}

		////At least k outlets, year 0+
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int k = data.params["station_coord"][j_bar][1];
				if (k > 1){
					model.add(x[t][j_bar] <= x[t][j_bar - 1]);
				}
			}
		}

		///Redundant, but prevents crash from unextracted variables
		model.add(theta_obj >= 0);
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
					model.add(x[t][j_bar] >= 0);
			}
			model.add(theta[t] >= 0);
		}


		//Set objective via proxy. Cut methods define theta_obj via theta variables
		model.add(IloMaximize(env, theta_obj));


		//Create warmstart solution (greedy)
		{IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		vector<vector<int>> sol;
		Greedy_Improved G;
		G.SetData(data);
		G.Solve(false);
		sol = G.Solution;
		ObjectiveValue = G.SolutionQuality;

		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				int val = 0;
				if (sol[t][j] >= k) {
					val = 1;
				}
				startVar.add(x[t][j_bar]);
				startVal.add(val);
			}

			startVar.add(theta[t]);
			IloNum val = 0;
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				if (sol[t][j] >= k) {
					cover += data.a[t].col(j_bar);
					val += data.Ps[t](j_bar);
				}
			}
			val += data.weights[t].dot((VectorXd)(cover.array() >= 1).matrix().cast<double>());
			startVal.add(val);
		}


		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
		}

		////Upper bound for theta (ensures bounded problem)

		for (int t = 0; t < T; t++) {
			IloNum ub = data.weights[t].sum() + data.Precovered[t] + data.Ps[t].sum();
			model.add(theta[t] <= ub);
		}

		//Objective
		if (verbose) { cplex.out() << "Adding objective \n"; }
		IloExpr obj(env);
		for (int t = 0; t < T; t++) { obj += theta[t]; }
		model.add(theta_obj == obj);
		obj.end();


		/////////////////////////////////////////////////////////////////////////

		//Link callback
		Callback_Improved cb(data, x, theta);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);


		bool skip_solve = false;
		////Trust region
		//if (use_trust) {
		//	IloConstraintArray trust(env);

		//	//Note: Trust distance definitions depend only on the warmstart solution, so never need to be updated
		//	IloExpr trust_expr(env);
		//	NumVar3D delta(env, T);
		//	//Num3D I_tilde(env, T);
		//	for (int t = 0; t < T; t++) {
		//		delta[t] = NumVar2D(env, M);
		//		//I_tilde[t] = Num2D(env, N);
		//		for (int j = 0; j < M; j++) {
		//			delta[t][j] = IloNumVarArray(env, Mj[j], 0.0, IloInfinity); //NumVar initialisation requires the bounds, or accessing it will crash!
		//			for (int k = 1; k < Mj[j]; k++) {
		//				IloExpr expr(env);
		//				expr -= delta[t][j][k];
		//				for (int t_bar = 0; t_bar < t; t_bar++) { expr += delta[t_bar][j][k]; } //We add the previous delta variables to account for outlet preservation constraints
		//				if (sol[t][j] >= k) { expr += 1 - x[t][j][k]; }
		//				else { expr += x[t][j][k]; }
		//				model.add(expr == 0);
		//				expr.end();
		//				trust_expr += delta[t][j][k];
		//			}
		//		}


		//	{
		//		IloConstraint cons = trust_expr <= trust_threshold;
		//		model.add(cons);
		//		trust.add(cons);
		//		//cons.end();
		//	}

		//	IloBool trust_solved = cplex.solve();
		//	if (trust_solved) {
		//		if (verbose) {
		//			cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
		//			cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
		//			cplex.out() << "Number of nodes: " << cplex.getNnodes() << endl;
		//			cplex.out() << "Trust region solve time: " << cplex.getTime() << endl;
		//		}
		//		ObjectiveValue = cplex.getObjValue();
		//		SolveTime = cplex.getTime();
		//		TotalTime = time(NULL) - start;
		//		OptimalityGap = cplex.getMIPRelativeGap();
		//		GetSolution(cplex, x);
		//		stats = cb.stats;
		//		stats["nNodes"] = cplex.getNnodes();
		//		stats["TrustSolveTime(DivideBy100)"] = (int)(100 * cplex.getTime());
		//	}
		//	//I_tilde.end();
		//	IloNum timelimit = 7200 - cplex.getTime();
		//	if (timelimit > EPS) { cplex.setParam(IloCplex::Param::TimeLimit, timelimit); }
		//	else { skip_solve = true; }
		//	//for (int i = 0; i < trust.getSize(); i ++) {
		//	//	model.remove(trust[i]);
		//	//}
		//	cplex.clearUserCuts();
		//	cplex.deleteMIPStarts(0, cplex.getNMIPStarts());
		//	{
		//		model.remove(trust);
		//		IloConstraint cons = trust_expr >= trust_threshold + 1;
		//		model.add(cons);
		//		trust.add(cons);
		//	}

		//	IloConstraint cons;
		//	IloExpr obj(env);
		//	switch (multicut)
		//	{

		//	case multicuts::SingleB1:
		//	case multicuts::SingleB2:
		//	case multicuts::SinglePO1:
		//		cons = theta[0][0] >= ObjectiveValue;
		//		model.add(cons);
		//		break;
		//	case multicuts::Multi1SB1:
		//	case multicuts::Multi1B1:
		//	case multicuts::Multi1B2:
		//	case multicuts::Multi1PO1:
		//		for (int t = 0; t < T; t++) { obj += theta[t][0]; }
		//		cons = obj >= ObjectiveValue;
		//		model.add(cons);
		//		break;
		//	case multicuts::Multi2B1:
		//	case multicuts::Multi2B2:
		//		for (int i = 0; i < N; i++) { obj += theta[0][i]; }
		//		cons = obj >= ObjectiveValue;
		//		model.add(cons);
		//		break;
		//	case multicuts::Multi3SB2:
		//	case multicuts::Multi3B1:
		//	case multicuts::Multi3B2:
		//	case multicuts::Multi3PO1:
		//		for (int t = 0; t < T; t++) { for (int i = 0; i < N; i++) { obj += theta[t][i]; } }
		//		cons = obj >= ObjectiveValue;
		//		model.add(cons);
		//		break;
		//	default:
		//		break;
		//	}
		//}

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

void Model_Improved::GetSolution(IloCplex& cplex, BoolVar2D& x)
{
	vector<vector<int>> temp;
	for (int t = 0; t < data.T; t++) {
		vector<int> Solution1;
		for (int j = 0; j < data.params["M"]; j++) { Solution1.push_back(0); }
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			if (cplex.getValue(x[t][j_bar]) > 1 - EPS) {
					Solution1[j] = k;
			}
		}
		temp.push_back(Solution1);
	}
	Solution = temp;
}







