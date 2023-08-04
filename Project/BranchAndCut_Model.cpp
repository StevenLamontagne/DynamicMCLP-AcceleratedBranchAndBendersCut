#include "BranchAndCut_Model.h"


void BranchAndCut_Model::Solve(json params)
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

	try {

	//Create model and set parameters
	IloModel model(env);
	IloCplex cplex(model);

	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::TimeLimit, 7200);
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
	

	BoolVar2D w(env, T);
	for (int t = 0; t < T; t++) {
		w[t] = IloBoolVarArray(env, P[t] + M_bar); //We add an extra M_bar to account for single-covered users, which the data groups separately 
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

	////Covering 
	for (int t = 0; t < T; t++) {
		//The data class is optimised for Benders methods, and so we must rebuild the original coverage matrix for use here
		IloArray<IloNumArray> coverage(env, P[t]);
		for (int p = 0; p < P[t]; p++) {
			coverage[p] = IloNumArray(env, M_bar);
		}
		for (int out = 0; out < data.a[t].outerSize(); out++) {
			for (SparseXXd::InnerIterator it(data.a[t], out); it; ++it) {
				int p = it.row();
				int j_bar = it.col();
				coverage[p][j_bar] = 1;
			}
		}


		for (int p = 0; p < P[t]; p++) {
			model.add(IloScalProd(coverage[p], x[t]) >= w[t][p]);
		}

		
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			model.add(x[t][j_bar] >= w[t][P[t] + j_bar]); //Single-covered users for each outlet
		}
	}


	//Objective
	IloExpr obj(env);
	for (int t = 0; t < T; t++) {
		obj += data.Precovered[t];
		for (int p = 0; p < P[t]; p++) {
			obj += data.weights[t](p) * w[t][p];
		}

		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			obj += data.Ps[t](j_bar) * w[t][P[t] + j_bar]; //Single-covered users for each outlet
		}
	}
	model.add(IloMaximize(env, obj));
	obj.end();


	cout << "Model creation time: " << time(NULL) - start << " seconds" << endl;



	//Solve and get results
	IloBool solved = cplex.solve();
	if (solved) {
		cplex.out() << "Solution status:" << cplex.getCplexStatus() << endl;
		cplex.out() << "Optimal value:" << cplex.getObjValue() << endl;


		stats["nNodes"] += (int)cplex.getNnodes();
		stats["CplexStatus"] = cplex.getStatus();
		stats["Solve time, MIP (x100)"] = (int)100 * cplex.getTime();
		stats["ObjectiveValue (x100)"] = (int)100 * cplex.getObjValue();
		stats["OptimalityGap (x100)"] = (int)100 * cplex.getMIPRelativeGap();

		GetSolution(cplex, x);
		ConvertSolution();
	}
	else {
		stats["ObjectiveValue (x100)"] = -1;
		stats["Solve time, MIP (x100)"] = -1;
		stats["OptimalityGap (x100)"] = -1;
		stats["nNodes"] = -1;
		stats["CplexStatus"] = 0;

		Solution = ArrayXXd::Constant(T, M_bar, 0.0);
		ConvertSolution();
	}


	}
	catch (IloException & e) {
		cout << "Exception: " << e << endl;
		stats["ObjectiveValue (x100)"] = -1;
		stats["Solve time, MIP (x100)"] = -1;
		stats["OptimalityGap (x100)"] = -1;
		stats["nNodes"] = -1;
		stats["CplexStatus"] = 0;

		Solution = ArrayXXd::Constant(data.T, data.M_bar, 0.0);
		ConvertSolution();
	}
	env.end();
}




