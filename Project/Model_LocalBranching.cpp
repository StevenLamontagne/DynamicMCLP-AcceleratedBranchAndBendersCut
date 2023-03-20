#include "Model_LocalBranching.h"

#define EPS 1.0e-6

ILOSTLBEGIN
//typedef IloArray<IloNumArray>    Num2D;
//typedef IloArray<IloNumVarArray>    NumVar2D;
//typedef IloArray<NumVar2D>    NumVar3D;
//typedef IloArray<IloBoolVarArray> BoolVar2D;
//typedef IloArray<BoolVar2D> BoolVar3D;


void Model_LocalBranching::SetData(const Data_Improved& newData)
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


void Model_LocalBranching::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);
	overlap_threshold = params.value("overlap_threshold", 0.1);
	use_trust = params.value("use_trust", false);
	bool use_trust_cuts = params.value("use_trust_cuts", false);
	trust_threshold = params.value("trust_threshold", 2.0);

	//string outfile = params.value("outfile", "");
	//ofstream fout(outfile);

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	Solution = ArrayXXd::Constant(T, M_bar, 0.0);
	//int M = data.params["M"];
	//vector<int> Mj;
	//for (int val : data.params["Mj"]) { Mj.push_back(val); }





	time_t start;
	time(&start);
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	VectorXd test;

	//if (outfile != "") {
	//	cplex.setOut(fout);
	//	cplex.setWarning(fout);
	//	cplex.setError(fout);
	//	cplex.out() << "Log sent to file" << endl;
	//}
	//else {
	//	cout << "No log file provided" << endl;
	//}
	try {
		//Set parameters
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200);

		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);

		cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 2);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::Param::Emphasis::MIP, 0);

		cplex.setParam(IloCplex::Param::MIP::Display, 1);
		cplex.setParam(IloCplex::Param::MIP::Interval, 10000);
		cplex.setParam(IloCplex::Param::Simplex::Display, 0);
		

		//Create variables
		BoolVar2D x(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = IloBoolVarArray(env, M_bar);
		}

		IloNumVarArray  theta(env, T, 0.0, IloInfinity);

		IloNumVar theta_obj(env);

		BoolVar2D trust(env, 1024);
		for (int i = 0; i < 1024; i++) {
			trust[i] = IloBoolVarArray(env, T);
		}


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
				if (k > 1) {
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

			for (int i = 0; i < 1024; i++) {
				model.add(trust[i][t] >= 0);
			}
		}


		//Create warmstart solution (greedy)
		IloNumVarArray startVar(env);
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
				Solution(t, j_bar) = val;
			}

			startVar.add(theta[t]);
			IloNum val = 0;
			VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				if (sol[t][j] >= k) {
					I_tilde += data.a[t].col(j_bar);
				}
			}
			val += data.weights[t].dot((VectorXd)(I_tilde.array() >= 1).matrix().cast<double>());
			startVal.add(val);
		}


		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();



		////Upper bound for theta (ensures bounded problem)
		for (int t = 0; t < T; t++) {
			model.add(theta[t] <= data.weights[t].sum());
		}


		//Objective
		if (verbose) { cplex.out() << "Adding objective \n"; }
		IloExpr obj(env);
		for (int t = 0; t < T; t++) {
			obj += theta[t] + data.Precovered[t];
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				obj += data.Ps[t](j_bar) * x[t][j_bar];
			}
		}
		model.add(theta_obj == obj);
		obj.end();

		//Set objective via proxy
		model.add(IloMaximize(env, theta_obj));


		//Link callback
		Callback_LocalBranching cb(data, model, x, theta, trust); 
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);


		////////////////////////////////////////////////////////////////////////////////////////////
		//Solve and get results
		bool solved = cplex.solve();
		if (verbose) {
			cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value: " << cb.GetSolutionObjective() << endl;
			//cplex.out() << "Number of nodes: " << cplex.getNnodes() << endl;
			cplex.out() << "\n" << endl;
		}

		for (pair<string, int> res : cb.stats) {
			string category = res.first;
			int value = res.second;
			stats[category] += value;
		}
		//stats["nNodes"] += (int)cplex.getNnodes();
		IloAlgorithm::Status status = cplex.getStatus();
		stats["CplexStatus"] = status;
		stats["SolveTime (x100)"] = (int)100 * cplex.getTime();
		stats["ObjectiveValue (x100)"] = (int)100 * cb.GetSolutionObjective();
		//stats["OptimalityGap (x100)"] = (int)10000 * cplex.getMIPRelativeGap();

		if (cb.LazyCutTimes.size() > 0) {
			stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
		}
		if (cb.UserCutTimes.size() > 0) {
			stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
		}

		Solution = cb.GetSolution();
		//stats["ObjectiveValue (x100)"] = 100 * cb.GetSolutionObjective();
	}
	catch (IloException & e) {
		int temp = 0;
		cplex.out() << "Exception: " << e << endl;
		ObjectiveValue = -1;
		SolveTime = -1;
		OptimalityGap = -1;
	}
	cplex.clearUserCuts();
	cplex.end();
	model.end();
	env.end();
	//fout.close();

}

void Model_LocalBranching::GetSolution(IloCplex& cplex, BoolVar2D& x)
{
	for (int t = 0; t < data.T; t++) {
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			Solution(t, j_bar) = cplex.getValue(x[t][j_bar]);
		}
	}
}