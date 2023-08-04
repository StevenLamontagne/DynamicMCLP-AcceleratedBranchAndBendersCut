#include "LocalBranching_Model.h"


void LocalBranching_Model::SetData(const Data& newData)
{
	data = newData;

}


void LocalBranching_Model::Solve(json params)
{
	//Set model parameters
	verbose = params.value("verbose", false);
	nTrust = params.value("nTrust", 1024);
	budgetType = params.value("budgetType", BUDGET_TYPE::Knapsack);



	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	
	//Initialise model and solution objects
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);
	Solution = ArrayXXd::Constant(T, M_bar, 0.0);
	ValueByYear = ArrayXd::Constant(T, 0.0);


	try {
		//Set parameters
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200);

		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);

		////Disable the presolve options that CPLEX bans with calbacks
		cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 1);

		cplex.setParam(IloCplex::Param::MIP::Limits::RepairTries, -1);
		cplex.setParam(IloCplex::Param::ParamDisplay, 1);
		cplex.setParam(IloCplex::Param::Read::WarningLimit, 0);
		cplex.setParam(IloCplex::Param::MIP::Display, 2);
		cplex.setParam(IloCplex::Param::MIP::Interval, 100);
		cplex.setParam(IloCplex::Param::Simplex::Display, 2);
		

		//Create variables
		x = BoolVar2D(env, T);
		for (int t = 0; t < T; t++) {
			x[t] = IloBoolVarArray(env, M_bar);
		}

		theta = IloNumVarArray(env, T, 0.0, IloInfinity);
		IloNumVar theta_obj(env);

		BoolVar2D trust(env, nTrust);
		for (int i = 0; i < nTrust; i++) {
			trust[i] = IloBoolVarArray(env, T);
		}


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
				if (k > 1) {
					model.add(x[t][j_bar] <= x[t][j_bar - 1]);
				}
			}
		}


		////Upper bound for theta (ensures bounded problem)
		for (int t = 0; t < T; t++) {
			model.add(theta[t] <= data.weights[t].sum());
		}


		////Ensures trust cut big-Ms work correctly
		////Without these constraints, it's unclear how CPLEX handles the trust cuts
		for (auto i = 0; i < trust.getSize(); i++) {
			model.add(IloSum(trust[i]) <= T - 1);
		}


		///Redundant, but prevents crash from unextracted variables
		model.add(theta_obj >= 0);
		for (int t = 0; t < T; t++) {
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				model.add(x[t][j_bar] >= 0);
			}
			model.add(theta[t] >= 0);

			for (auto i = 0; i < trust.getSize(); i++) {
				model.add(trust[i][t] >= 0);
			}
		}

		//Objective
		IloExpr obj(env);
		for (int t = 0; t < T; t++) {
			obj += theta[t] + data.Precovered[t];
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				obj += data.Ps[t](j_bar) * x[t][j_bar];
			}
		}
		model.add(theta_obj == obj);
		obj.end();

		model.add(IloMaximize(env, theta_obj)); //Set objective via proxy


		//Calculate solution of LP relaxation and add associated fractional Benders cuts to model
		{
			vector<IloConversion> conv_x;
			for (int t = 0; t < T; t++) {
				IloConversion conv(env, x[t], IloNumVar::Type::Float);
				model.add(conv);
				conv_x.push_back(conv);
			}

			MultiCutBenders_Callback cb(data, x, theta);
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation;
			cplex.use(&cb, contextmask);
			cplex.solve();
			stats["Solve time, LP (x100)"] = (int)100 * cplex.getTime();
			stats["Objective value, LP (x100)"] = (int)100 * cplex.getObjValue();
			cout << "Objective value, LP: " << (double) stats["Objective value, LP (x100)"]/ (double) 100.00 << endl;

			for (auto conv : conv_x) {
				model.remove(conv);
			}
			cplex.deleteMIPStarts(0);

			cplex.use(NULL, 0);


		}

		
		//Create warmstart solution (greedy)
		{
			vector<vector<int>> sol;
			Greedy G;
			G.SetData(data);
			G.Solve(false, budgetType);
			sol = G.Solution;
			if (G.SolutionQuality > ObjectiveValue) {
				ObjectiveValue = G.SolutionQuality;
				for (int t = 0; t < T; t++) {
					for (int j_bar = 0; j_bar < M_bar; j_bar++) {
						int j = data.params["station_coord"][j_bar][0];
						int k = data.params["station_coord"][j_bar][1];
						int val = 0;
						if (sol[t][j] >= k) {
							val = 1;
						}
						Solution(t, j_bar) = val;
					}
				}

			}
		}


		//Link local branching callback
		//Solving takes place here too for scoping reasons


		//To use the SepB method, uncomment the following block
		//////////////////////////////////////////////////////////////////////////////////

		//if (T == 2) {
		//	LocalBranching_TwoYearCallback cb(data, model, x, theta, theta_obj, trust);
		//	CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
		//		| IloCplex::Callback::Context::Id::Relaxation
		//		| IloCplex::Callback::Context::Id::Branching;
		//	cplex.use(&cb, contextmask);
		//	cb.SetSolution(Solution, ObjectiveValue);

		//	//Solve and get results
		//	bool solved = cplex.solve();
		//	cplex.out() << endl;
		//	if (verbose) {
		//		cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
		//		cplex.out() << "Optimal value: " << cb.GetSolutionObjective() << endl;
		//		cplex.out() << "\n" << endl;
		//	}

		//	for (pair<string, int> res : cb.stats) {
		//		string category = res.first;
		//		int value = res.second;
		//		stats[category] = value;
		//	}

		//	stats["nNodes"] += (int)cplex.getNnodes();
		//	stats["CplexStatus"] = (int)cplex.getStatus();
		//	stats["Solve time, MIP (x100)"] = (int)100 * cplex.getTime();
		//	stats["ObjectiveValue (x100)"] = (int)100 * cb.GetSolutionObjective();
		//	stats["OptimalityGap (x100)"] = (int)10000 * cplex.getMIPRelativeGap();

		//	if (cb.LazyCutTimes.size() > 0) {
		//		stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
		//	}
		//	if (cb.UserCutTimes.size() > 0) {
		//		stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
		//	}

		//	Solution = cb.GetSolution();
		//}
		//else{
			LocalBranching_Callback cb(data, model, x, theta, theta_obj, trust);
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
				| IloCplex::Callback::Context::Id::Relaxation;
			cplex.use(&cb, contextmask);
			cb.SetSolution(Solution, ObjectiveValue);

			//Solve and get results
			bool solved = cplex.solve();
			cplex.out() << endl;
			if (verbose) {
				cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
				cplex.out() << "Optimal value: " << cb.GetSolutionObjective() << endl;
				cplex.out() << "\n" << endl;
			}

			for (pair<string, int> res : cb.stats) {
				string category = res.first;
				int value = res.second;
				stats[category] = value;
			}
			stats["nNodes"] += (int)cplex.getNnodes();
			stats["CplexStatus"] = (int)cplex.getStatus();
			stats["Solve time, MIP (x100)"] = (int)100 * cplex.getTime();
			stats["ObjectiveValue (x100)"] = (int)100 * cb.GetSolutionObjective();
			stats["OptimalityGap (x100)"] = (int)10000 * cplex.getMIPRelativeGap();

			if (cb.LazyCutTimes.size() > 0) {
				stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
			}
			if (cb.UserCutTimes.size() > 0) {
				stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
			}

			Solution = cb.GetSolution();
		
		//} 

		////////////////////////////////////////////////////////////////////////////////////////////
		//End of the block for SepB method



		ConvertSolution();

	}
	catch (IloException & e) {
		int temp = 0;
		cplex.out() << "Exception: " << e << endl;
		ObjectiveValue = -1;
	}

	env.end();
}


