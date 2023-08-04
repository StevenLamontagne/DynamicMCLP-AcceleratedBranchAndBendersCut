#include "MultiCutBenders_Model.h"


void MultiCutBenders_Model::SetData(const Data& newData)
{
	data = newData;

}


void MultiCutBenders_Model::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);
	budgetType = params.value("budgetType", BUDGET_TYPE::Knapsack);


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


		IloNumVarArray  theta(env, T, 0.0, IloInfinity);

		IloNumVar theta_obj(env);


		//Constraints
		switch (budgetType)
		{
		case BUDGET_TYPE::Knapsack:
		{
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
			break;


		}
		case BUDGET_TYPE::OutletCount:
		{
			////Max outlets installed per time period, year 0 
			IloExpr outlets0(env);
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				outlets0 += x[0][j_bar] - (int)data.params["x0"][j][k];
			}
			model.add(outlets0 <= (int)data.params["Stations_maxNewOutletsPerTimePeriod"]);
			outlets0.end();

			IloExpr stations0(env);
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				if (k == 1) { stations0 += x[0][j_bar] - (int)data.params["x0"][j][k]; }
			}
			model.add(stations0 <= (int)data.params["Stations_maxNewStationsPerTimePeriod"]);
			stations0.end();

			////Max outlets installed per time period, year 1+ 
			for (int t = 1; t < T; t++) {
				IloExpr outlets(env);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) {
					outlets += (x[t][j_bar] - x[t - 1][j_bar]);
				}
				model.add(outlets <= (int)data.params["Stations_maxNewOutletsPerTimePeriod"]);
				outlets.end();

				IloExpr stations(env);
				for (int j_bar = 0; j_bar < M_bar; j_bar++) {
					int k = data.params["station_coord"][j_bar][1];
					if (k == 1) { stations += x[t][j_bar] - x[t - 1][j_bar]; }
				}
				model.add(stations <= (int)data.params["Stations_maxNewStationsPerTimePeriod"]);
				stations.end();

			}
			break;
		}
		default:
		{
			cout << "Unrecognised budget type." << endl;
			throw;
			break;
		}
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



		////////////////////////////////////////////////////////////////////////////////////////////
		//Create warmstart solution (greedy)
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		vector<vector<int>> sol;
		Greedy G;
		G.SetData(data);
		G.Solve(false, budgetType);
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
		//////////////////////////////////////////////////////////////////////////////////////////////


		////Upper bound for theta (ensures bounded problem)

		for (int t = 0; t < T; t++) {
			IloNum ub = data.weights[t].sum() + data.Precovered[t] + data.Ps[t].sum();
			model.add(theta[t] <= ub);
		}

		////////////////////////////////////////////////////////////////////////////////////////////

		//Objective
		if (verbose) { cplex.out() << "Adding objective \n"; }
		//Objective
		IloExpr obj(env);
		for (int t = 0; t < T; t++) {
			obj += theta[t] + data.Precovered[t];
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				obj += data.Ps[t](j_bar) * x[t][j_bar];
			}
		}
		model.add(IloMaximize(env, obj));
		obj.end();


		///////////////////////////////////////////////////////////////////////////

		////Link callback
		MultiCutBenders_Callback cb(data, x, theta); 
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);
		


		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t];
			lhs += data.Precovered[t];

			for (int j_bar = 0; j_bar < M_bar;j_bar++) {
				lhs += (data.CutCoeffs[t].col(j_bar).sum() + data.Ps[t](j_bar)) * x[t][j_bar];
			}
			model.add(lhs >= 0);
			lhs.end();
		}

		
		//Solve and get results
		cplex.solve();
		IloAlgorithm::Status status = cplex.getStatus();
		if (verbose) {
			cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
			cplex.out() << "\n" << endl;
		}

		for (pair<string, int> res : cb.stats) {
			string category = res.first;
			int value = res.second;
			stats[category] += value;
		}
		stats["nNodes"] += (int) cplex.getNnodes();
		stats["CplexStatus"] = status;
		stats["Solve time, MIP (x100)"] = (int) 100* cplex.getTime();
		stats["ObjectiveValue (x100)"] = (int) 100* cplex.getObjValue();
		stats["OptimalityGap (x100)"] = (int) 10000 * cplex.getMIPRelativeGap();
		
		if (cb.LazyCutTimes.size() > 0) {
			stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
		}
		if (cb.UserCutTimes.size() > 0) {
			stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
		}
		

	}
	catch (IloException & e) {
		cplex.out() << "Exception: " << e << endl;
		ObjectiveValue = -1;		
		stats["nNodes"] += -1;
		stats["CplexStatus"] = 0;
		stats["Solve time, MIP (x100)"] = -1;
		stats["ObjectiveValue (x100)"] = -1;
		stats["OptimalityGap (x100)"] = -1;

	}
	cplex.clearUserCuts();
	cplex.end();
	model.end();
	env.end();
	//fout.close();

}







