#include "SingleCutBenders_Model.h"

void SingleCutBenders_Model::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;


	Solution = ArrayXXd::Constant(T, M_bar, 0.0);


	time_t start;
	time(&start);


	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);
	try {
		//Set parameters
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200);
		cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 2);
		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);
		cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);
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

		IloNumVar theta(env);


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

		////Upper bound for theta (ensures bounded problem)
		IloNum ub = 0;
		for (int t = 0; t < T; t++) { ub += data.weights[t].sum() + data.Ps[t].sum(); }
		model.add(theta <= ub);


		///Redundant, but prevents crash from unextracted variables
		model.add(theta >= 0);
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				model.add(x[t][j_bar] >= 0);
			}
		}

		//Objective
		cout << "Adding objective \n";
		IloExpr obj(env);
		obj += theta;
		for (int t = 0; t < T; t++ ){ obj += data.Precovered[t]; }
		model.add(IloMaximize(env, obj));
		obj.end();


		//Link callback
		SingleCutBenders_Callback cb(data, x, theta);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);




		cout << "Model creation time: " << time(NULL) - start << " seconds" << endl;

		//Solve and get results
		IloBool solved = cplex.solve();
		if (solved) {
			cplex.out() << "Solution status:" << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value:" << cplex.getObjValue() << endl;

			for (pair<string, int> res : cb.stats) {
				string category = res.first;
				int value = res.second;
				stats[category] = value;
			}
			stats["nNodes"] += (int)cplex.getNnodes();
			stats["CplexStatus"] = (int)cplex.getStatus();
			stats["Solve time, MIP (x100)"] = (int)100 * cplex.getTime();
			stats["ObjectiveValue (x100)"] = (int)100 * cplex.getObjValue();
			stats["OptimalityGap (x100)"] = (int)10000 * cplex.getMIPRelativeGap();

			if (cb.LazyCutTimes.size() > 0) {
				stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
			}
			if (cb.UserCutTimes.size() > 0) {
				stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
			}

			GetSolution(cplex, x);
			ConvertSolution();


		}
		else {
			ObjectiveValue = -1;
			stats["nNodes"] = -1;
			stats["CplexStatus"] = 0;
			stats["Solve time, MIP (x100)"] = -1;
			stats["ObjectiveValue (x100)"] =  -1;
			stats["OptimalityGap (x100)"] = -1;
		}
	}
catch (IloException & e) {
	cout << "Exception: " << e << endl;
	ObjectiveValue = -1;
	stats["nNodes"] = -1;
	stats["CplexStatus"] = 0;
	stats["Solve time, MIP (x100)"] = -1;
	stats["ObjectiveValue (x100)"] = -1;
	stats["OptimalityGap (x100)"] = -1;
}
env.end();
}

