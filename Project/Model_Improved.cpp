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
	//overlap_threshold = params.value("overlap_threshold", 0.1);
	//use_trust = params.value("use_trust", false);
	//bool use_trust_cuts = params.value("use_trust_cuts", false);
	//trust_threshold = params.value("trust_threshold", 2.0);

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
		//Create model and set parameters
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, 7200); //REMEMBER TO RESET TO 7200!!!!!!!!!!!!!!

		//cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
		//cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 2);

		cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		cplex.setParam(IloCplex::Param::WorkMem, 100000);
		cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);

		cplex.setParam(IloCplex::Param::MIP::Display, 1);
		cplex.setParam(IloCplex::Param::MIP::Interval, 10000);
		cplex.setParam(IloCplex::Param::Simplex::Display, 0);
	
		//cplex.setParam(IloCplex::Param::Network::Display, 0);
		//cplex.setParam(IloCplex::Param::Sifting::Display, 0);		
		//cplex.setParam(IloCplex::Param::Barrier::Display, 0);
		//cplex.setParam(IloCplex::Param::Conflict::Display, 0);
		//cplex.setParam(IloCplex::Param::Network::Display, 0);	
		//cplex.setParam(IloCplex::Param::Tune::Display, 0);


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

		////////////////////////////////////////////////////////////////////////////////////////////
		//Set objective via proxy. Cut methods define theta_obj via theta variables
		model.add(IloMaximize(env, theta_obj));

		////////////////////////////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////////////////////////////
		cout << "Adding trust constraints" << endl;
		

		////Trust constraints for greedy solution (distance 2)
		//for (int t = 0; t < T; t++) {
		//	IloExpr trust_expr(env);
		//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//		if (Solution(t, j_bar) > 1 - EPS) { trust_expr += 1 - x[t][j_bar]; }
		//		else { trust_expr += x[t][j_bar]; }
		//	}
		//	model.add(trust_expr <= 2);
		//}

		//Trust constraints for greedy solution (distance 4)
		IloExpr taboo(env);
		for (int t = 0; t < T; t++) {
			IloExpr trust_expr(env);

			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				if (Solution(t, j_bar) > 1 - EPS) { trust_expr += 1 - x[t][j_bar]; }
				else { trust_expr += x[t][j_bar]; }
			}
			model.add(trust_expr <= 4);
			taboo += trust_expr;
		}
		model.add(taboo >= 1);
		taboo.end();

		////////////////////////////////////////////////////////////////////////////////////////////

		////Two-opt cut generation
		//for (int t = 0; t < T; t++) {
		//	vector<int> always;
		//	vector<int> active;
		//	vector<int> inactiveD1;
		//	vector<int> inactiveD2;
		//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//		int k = data.params["station_coord"][j_bar][1];
		//		if (Solution(t, j_bar) > 1 - EPS) {
		//			int j = data.params["station_coord"][j_bar][0];
		//			if ((k == ((int) data.params["Mj"][j] - 1)) || (Solution(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
		//			else { 
		//				always.push_back(j_bar); 
		//				model.add(x[t][j_bar] == Solution(t, j_bar));
		//			}
		//		}
		//		else if ((k == 1) || (Solution(t, j_bar - 1) > 1 - EPS)) { inactiveD1.push_back(j_bar); }
		//		else if ((k == 2) || (Solution(t, j_bar - 2) > 1 - EPS)) { inactiveD2.push_back(j_bar); }
		//		else{ model.add(x[t][j_bar] == Solution(t, j_bar)); }

		//	}

		//	//Initialise coverage
		//	VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
		//	for (int j_bar : always) { I_tilde += data.a[t].col(j_bar); }
		//	for (int j_bar : active) { I_tilde += data.a[t].col(j_bar); }

		//	//Cut for greedy solution + distance 1
		//	{
		//		IloExpr lhs(env);
		//		lhs -= theta[t];
		//		lhs += (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		//		VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
		//		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
		//			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		//		}
		//		model.add(lhs >= 0);
		//	}

		//	//Cuts for adding two outlets (redundant for greedy solution)
		//	for (int j_bar : inactiveD1) {

		//		IloExpr lhs(env);
		//		lhs -= theta[t];

		//		VectorXd mod_I_tilde = I_tilde + data.a[t].col(j_bar);
		//		lhs += (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		//		VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
		//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		//		for (int j_hat : inactiveD2) {
		//			lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
		//		}
		//		model.add(lhs >= 0);
		//	}

		//	//Check coverage when swapping out active stations with different ones
		//	for (int j_bar : active) {

		//		IloExpr lhs(env);
		//		lhs -= theta[t];

		//		VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
		//		lhs += (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
		//		VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
		//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		//		for (int j_hat : inactiveD1) {
		//			lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
		//		}
		//		model.add(lhs >= 0);
		//	}
		//}

		////////////////////////////////////////////////////////////////////////////////////////////

		////Upper bound for theta (ensures bounded problem)

		for (int t = 0; t < T; t++) {
			model.add(theta[t] <= data.weights[t].sum());
		}

		///////////////////////////////////////////////////////////////////////////////////////////


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


		///////////////////////////////////////////////////////////////////////////

		////Link callback
		Callback_Improved cb(data, x, theta); //Make sure only base model is included when running this
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::Relaxation;
		cplex.use(&cb, contextmask);
		

		//Solve and get results
		cplex.solve();
		//cplex.out() << "Max Primal residual: " << cplex.getQuality(IloCplex::Quality::MaxPrimalResidual) << endl;
		//cplex.out() << "Max Primal infeasibility: " << cplex.getQuality(IloCplex::Quality::MaxPrimalInfeas) << endl;
		//cplex.out() << "Kappa stable: " << cplex.getQuality(IloCplex::Quality::KappaStable) << endl;
		//cplex.out() << "Kappa suspicious: " << cplex.getQuality(IloCplex::Quality::KappaSuspicious) << endl;
		//cplex.out() << "Kappa unstable: " << cplex.getQuality(IloCplex::Quality::KappaUnstable) << endl;
		//cplex.out() << "Kappa ill-posed: " << cplex.getQuality(IloCplex::Quality::KappaIllposed) << endl;
		//cplex.out() << "Kappa attention: " << cplex.getQuality(IloCplex::Quality::KappaAttention) << endl;
		if (verbose) {
			cplex.out() << "Solution status: " << cplex.getCplexStatus() << endl;
			cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
			cplex.out() << "Number of nodes: " << cplex.getNnodes() << endl;
			cplex.out() << "\n" << endl;
		}

		for (pair<string, int> res : cb.stats) {
			string category = res.first;
			int value = res.second;
			stats[category] += value;
		}
		stats["nNodes"] += (int) cplex.getNnodes();
		IloAlgorithm::Status status = cplex.getStatus();
		stats["CplexStatus"] = status;
		stats["SolveTime (x100)"] = (int)100 * cplex.getTime();
		stats["ObjectiveValue (x100)"] = (int) 100* cplex.getObjValue();
		stats["OptimalityGap (x100)"] = (int) 10000* cplex.getMIPRelativeGap();
		
		if (cb.LazyCutTimes.size() > 0) {
			stats["LazyCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.LazyCutTimes.data(), cb.LazyCutTimes.size()).mean();
		}
		if (cb.UserCutTimes.size() > 0) {
			stats["UserCutTime (x1000)"] = 1000 * Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cb.UserCutTimes.data(), cb.UserCutTimes.size()).mean();
		}
		
		ObjectiveValue = cplex.getObjValue();
		GetSolution(cplex, x);
		OptimalityGap = cplex.getMIPRelativeGap();
		SolveTime = cplex.getTime();
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
	//fout.close();

}

void Model_Improved::GetSolution(IloCplex& cplex, BoolVar2D& x)
{
	for (int t = 0; t < data.T; t++) {
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			Solution(t, j_bar) = cplex.getValue(x[t][j_bar]);
		}
	}
}

vector<vector<int>> Model_Improved::ConvertSolution()
{
	vector<vector<int>> sol;
	for (int t = 0; t < data.T; t++) {
		vector<int> temp;
		for (int j = 0; j < data.params["M"]; j++) { temp.push_back(0); }
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			int j = data.params["station_coord"][j_bar][0];
			int k = data.params["station_coord"][j_bar][1];
			if (Solution(t, j_bar) > 1 - EPS) {
				temp[j] = k;
			}
		}
		sol.push_back(temp);
	}
	return sol;
}


//void Model_Improved::GetSolution(IloCplex& cplex, BoolVar2D& x)
//{
//	vector<vector<int>> temp;
//	for (int t = 0; t < data.T; t++) {
//		vector<int> Solution1;
//		for (int j = 0; j < data.params["M"]; j++) { Solution1.push_back(0); }
//		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
//			int j = data.params["station_coord"][j_bar][0];
//			int k = data.params["station_coord"][j_bar][1];
//			if (cplex.getValue(x[t][j_bar]) > 1 - EPS) {
//				Solution1[j] = k;
//			}
//		}
//		temp.push_back(Solution1);
//	}
//	Solution = temp;
//}





