#include "MulticutCallback.h"
#define EPS 1.0e-6



void MulticutCallback::invoke(const IloCplex::Callback::Context& context)
{
	//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	IloEnv env = context.getEnv();
	IloArray<IloNumArray> x_tilde(env, T);
	

	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
	{
		//Get current solution value, x_tilde
		IloArray<IloArray<IloNumArray>> x_tilde(env, T);
		for (int t = 0; t < T; t++) {
			x_tilde[t] = IloArray<IloNumArray>(env, M);
			for (int j = 0; j < M; j++) {
				x_tilde[t][j] = IloNumArray(env, Mj[j]);
				context.getCandidatePoint(x[t][j], x_tilde[t][j]);
			}
		}

		LazyCutCallback(context, x_tilde);
		//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
		if ((heuristic == useHeuristic::WarmstartAndPostSimple || heuristic == useHeuristic::WarmstartAndPostGreedy) && (context.getIncumbentObjective() != incumbent)){
			incumbent = (double) context.getIncumbentObjective();
			AddPostHeuristic(context, x_tilde);
		}
		//cout << "Lazy cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		//Get current solution value, x_tilde
		IloArray<IloArray<IloNumArray>> x_tilde(env, T);
		for (int t = 0; t < T; t++) {
			x_tilde[t] = IloArray<IloNumArray>(env, M);
			for (int j = 0; j < M; j++) {
				x_tilde[t][j] = IloNumArray(env, Mj[j]);
				context.getRelaxationPoint(x[t][j], x_tilde[t][j]);
			}
		}


		UserCutCallback(context, x_tilde);
		if ((heuristic == useHeuristic::WarmstartAndPostSimple || heuristic == useHeuristic::WarmstartAndPostGreedy) && (context.getIncumbentObjective() != incumbent)) {
			incumbent = (double) context.getIncumbentObjective();
			AddPostHeuristic(context, x_tilde);
		}
		//cout << "User cut time: " << chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count() << endl;

		break;
	}
	}
}


void MulticutCallback::LazyCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde)
{
	//time_t start;
	//time(&start);

	IloEnv env = context.getEnv();

	//Calculate I_tilde
	Num3D I_tilde(env, T);
	for (int t = 0; t < T; t++) {
		I_tilde[t] = IloArray<IloArray<IloNum>>(env, N);
		for (int i = 0; i < N; i++) {
			I_tilde[t][i] = IloArray<IloNum>(env, R[i]);
			for (int r = 0; r < R[i]; r++) {
				IloInt val = 0;
				for (auto cover : data.cover[t][i][r]) {
					int j = cover.first;
					int k0 = cover.second;
					for (int k = k0; k < Mj[j]; k++) {
						if (x_tilde[t][j][k] > 1-EPS) { val += 1; }
					}		
				}
				I_tilde[t][i][r] = val;
			}
		}
	}

	
	//if (cut_type == multicuts::Multi1SB2) { AddStrengthenedCuts_Lazy(context, x_tilde, I_tilde); }
	//else { 
	AddCuts(context, I_tilde); 
	//}
	I_tilde.end();
}

void MulticutCallback::UserCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde)
{
	//time_t start;
	//time(&start);

	IloEnv env = context.getEnv();


	//Calculate I_tilde
	Num3D I_tilde(env, T);
	for (int t = 0; t < T; t++) {
		I_tilde[t] = IloArray<IloArray<IloNum>>(env, N);
		for (int i = 0; i < N; i++) {
			I_tilde[t][i] = IloArray<IloNum>(env, R[i]);
			for (int r = 0; r < R[i]; r++) {
				IloNum val = 0;
				for (int j = 0; j < M; j++) {
					for (int k = 1; k < Mj[j]; k++) {
						if ((x_tilde[t][j][k] > EPS) && data.a[t][i][r][j][k]) {
							val += x_tilde[t][j][k];
						}
					}
				}

				I_tilde[t][i][r] = val;
			}
		}
	}


	if (cut_type == multicuts::Multi1SB2) { AddStrengthenedCuts_User(context, x_tilde, I_tilde); }
	else { AddCuts(context, I_tilde); }

	I_tilde.end();

}


pair<int, int> MulticutCallback::argmax(map<pair<int, int>, double> vec)
{
		pair<int, int> maxKey;
		double maxVal = 0.0;
		for (auto f : vec) {
			double val = f.second;
			if (val > maxVal) { 
				maxVal = val;
				maxKey = f.first;
			}
		}
		return maxKey;
}

vector<map<int, vector<int>>> MulticutCallback::GetFractional(const IloArray<IloArray<IloNumArray>>& x_tilde)
{
	vector<map<int, vector<int>>> fracList;
	for (int t = 0; t < data.T; t++) {
		map<int, vector<int>> fracMap;
		for (int j = 0; j < data.M; j++) {
			vector<int> frac;
			for (int k = 1; k < data.Mj[j]; k++) {
				if ((x_tilde[t][j][k] > EPS) && (x_tilde[t][j][k] < 1 - EPS)) { frac.push_back(k); }
			}
			if (frac.size() > 0) { fracMap[j] = frac; }
		}
		fracList.push_back(fracMap);
	}
	return fracList;
}


bool MulticutCallback::SimpleRepair(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage)
{
	auto frac = GetFractional(x_tilde);

	IloArray<IloArray<IloNumArray>> newSolution(env, data.T);
	for (int t = 0; t < data.T; t++) {
		newSolution[t] = IloArray<IloNumArray>(env, data.M);
		for (int j = 0; j < data.M; j++) {
			newSolution[t][j] = IloNumArray(env, data.Mj[j]);
			for (int k = 0; k < data.Mj[j]; k++) {
				if (x_tilde[t][j][k] > 1 - EPS) { newSolution[t][j][k] = 1.0; }
				else { newSolution[t][j][k] = 0.0; }
			}
		}
	}

	for (int t = 0; t < data.T; t++) {
		//Can skip year if no fractional outlets
		if (frac[t].size() == 0) { continue; }

		//Refund budget for fractional solutions and set values to 0
		//cout << "Initial budget: " << Budget[t] << endl;
		for (auto frac_station : frac[t]) {
			int j = frac_station.first;
			double y_old = 0.0;
			//if (t == 0) { cout << "Station " << j << ", initial: " << data.params["x0"][j] << endl; }
			//else { cout << "Station " << j << ", year " << t - 1 << ": " << x_tilde[t - 1][j] << endl; }
			//cout << "Station " << j << ", year " << t << ": " << x_tilde[t][j] << endl;

			if (t == 0) {
				y_old = (double)data.params["y0"][j];
			}
			else {
				for (int k = 1; k < data.Mj[j]; k++) {
					y_old += min(x_tilde[t - 1][j][k], 1.0);
				}
			}
			for (int k = 1; k < Mj[j]; k++) {
				double previous = 0.0;
				if ((t == 0) && (k == (int)data.params["x0"][j])) { previous = 1; }
				else if (t == 0) { previous = 0; }
				else { previous = x_tilde[t - 1][j][k]; }

				double delta = x_tilde[t][j][k] - previous;
				Budget[t] += delta * k * (double)data.params["c"][t][j];
			}
			double delta = -y_old;
			for (int k = 1; k < data.Mj[j]; k++) {
				delta += x_tilde[t][j][k];
			}
			if (delta > EPS) { Budget[t] += delta * (double)data.params["f"][t][j]; }
			//cout << "Adjusted budget: " << Budget[t] << endl;

		}

		//Simple repairing
		for (auto frac_station : frac[t]) {
			int j = frac_station.first;
			int k0 = frac_station.second[0];
			int previous = 0;
			if (t == 0) { previous = (int)data.params["x0"][j]; }
			else {
				for (int k = 1; k < Mj[j]; k++) {
					if (newSolution[t - 1][j][k] > 1 - EPS) { previous = k; }
				}
			}
			k0 = max(k0, previous);
			double cost = (double)data.params["c"][t][j] * (k0 - previous);
			if (previous == 0) { cost += (double)data.params["f"][t][j]; }

			if (k0 == previous) {newSolution[t][j][k0] = 1;}
			if ((k0 > previous) && (Budget[t] >= cost)) {
				newSolution[t][j][k0] = 1;
				if (previous == 0) { Budget[t] -= cost; }
			}
			else { newSolution[t][j][previous] = 1; }
		}

		if (Budget[t] < - EPS) {
			if (frac[t].size() == 0) { continue; }
			else { newSolution.end();  return false; }
		}
		
	}

	x_tilde = newSolution.copy();
	newSolution.end();
	return true;
}

bool MulticutCallback::GreedyRepair(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage)
{
	auto frac = GetFractional(x_tilde);

	IloArray<IloArray<IloNumArray>> newSolution(env, data.T);
	for (int t = 0; t < data.T; t++) {
		newSolution[t] = IloArray<IloNumArray>(env, data.M);
		for (int j = 0; j < data.M; j++) {
			newSolution[t][j] = IloNumArray(env, data.Mj[j]);
			for (int k = 0; k < data.Mj[j]; k++) {
				if (x_tilde[t][j][k] > 1 - EPS) { newSolution[t][j][k] = 1.0; }
				else {	newSolution[t][j][k] = 0.0;	}
			}
		}
	}
	
	for (int t = 0; t < data.T; t++) {
		//Can skip year if no fractional outlets
		if (frac[t].size() == 0) { continue; }

		//Refund budget for fractional solutions and set values to 0
		//cout << "Initial budget: " << Budget[t] << endl;
		for (auto frac_station : frac[t]) {
			int j = frac_station.first;
			double y_old = 0.0;
			//if (t == 0) { cout << "Station " << j << ", initial: " << data.params["x0"][j] << endl; }
			//else { cout << "Station " << j << ", year " << t -1<< ": " << x_tilde[t-1][j] << endl; }
			//cout << "Station " << j << ", year " << t << ": " << x_tilde[t][j] << endl;

			if (t == 0) { 
				y_old = (double) data.params["y0"][j]; 
			}
			else {
				for (int k = 1; k < data.Mj[j]; k++) {
					y_old += min(x_tilde[t - 1][j][k], 1.0);
			}
			}
			for (int k = 1; k < Mj[j];k++) {
				double previous = 0.0;
				if ((t == 0) && (k == (int) data.params["x0"][j])){ previous = 1; }
				else if (t == 0) { previous = 0; }
				else { previous = x_tilde[t - 1][j][k]; }

				double delta = x_tilde[t][j][k] - previous;
				Budget[t] += delta * k *(double)data.params["c"][t][j];
			}
			double delta = -y_old;
			for (int k = 1; k < data.Mj[j]; k++) {
				delta += x_tilde[t][j][k];
			}
			if (delta > EPS) { Budget[t] += delta * (double)data.params["f"][t][j]; }
			//cout << "Adjusted budget: " << Budget[t] << endl;

		}

		//Greedy repairing
		if (Budget[t] <= EPS) { 
			if (frac[t].size() == 0) { continue; }
			else { return false; }
		}

		map<pair<int,int>, double> StationTotals;
		map<pair<int, int>, double> Costs;
		bool foundImprovement;
		do
		{
			foundImprovement = false;
			//Reset totals
			for (auto frac_station : frac[t]) {
				int j = frac_station.first;
				for (int k : frac_station.second) {
					StationTotals[make_pair(j, k)] = 0.0;
					Costs[make_pair(j, k)] = 0.0;
				}
			}

			//Check for coverage from new stations
			for (int i = 0; i < data.N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
				for (int r = 0; r < data.R[i]; r++) {
					if (coverage[t][i][r]) { continue; }
					for (auto frac_station : frac[t]) {
						int j = frac_station.first;
						for (int k : frac_station.second) {
							if (data.a[t][i][r][j][k]) {
								StationTotals[make_pair(j, k)] +=  weight;
							}
						}
					}
				}
			}

			//Reset coverage of stations to 0 if infeasible
			for (auto frac_station : frac[t]) {
				int j = frac_station.first;
				double previous_y = 0.0;
				double previous_x = 0.0;
				if (t == 0) { previous_y = (double)data.params["x0"][j]; }
				else {
					for (int k = 1; k < Mj[j]; k++) {
						previous_y += newSolution[t - 1][j][k];
						previous_x += k* newSolution[t - 1][j][k];
					}
				}
				for (int k : frac_station.second) {
					//Number of outlets decreases
					if (k < previous_x) { StationTotals[make_pair(j, k)] = 0.0;  continue; }
					//Budget insufficient for new outlet
					double cost = (1 - previous_y) * (double)data.params["f"][t][j] + k * (double)data.params["c"][t][j];
					if (cost > Budget[t] + EPS) {
						StationTotals[make_pair(j, k)] = 0.0;
					}
					else {
						StationTotals[make_pair(j, k)] = StationTotals[make_pair(j, k)] / cost;
						Costs[make_pair(j, k)] = cost;
					}
				}
			}

			pair<int, int> jBar_star = argmax(StationTotals);
			double zBar_star = StationTotals[jBar_star];
			if (zBar_star > 0) {
				//Update budget
				Budget[t] = max(Budget[t] - Costs[jBar_star], 0.0);
				int j = jBar_star.first;
				int k = jBar_star.second;
				for (int tBar = t; tBar < data.T; tBar++) {
					//Update solution for current (and all future) years
					newSolution[tBar][j][k] = 1;

					//Update coverage of triplets
					for (int i = 0; i < data.N; i++) {
						for (int r = 0; r < data.R[i]; r++) {
							if ((coverage[tBar][i][r]) && (data.a[tBar][i][r][j][k])) {
								coverage[tBar][i][r] = 1;
							}
						}
					}
				}
				foundImprovement = 1;
				frac[t].erase(j);
				StationTotals.clear();
				Costs.clear();
				if (frac[t].size() == 0) { break; }
			}
			else {
				if (frac[t].size() > 0) { return false; }
			}
		} while (foundImprovement);
	}


	x_tilde = newSolution.copy();
	newSolution.end();
	return true;
}

bool MulticutCallback::GreedyFill(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage)
{	//Initialise vectors for holding information for each station
	map<pair<int, int>, double> StationTotals;
	map<pair<int, int>, double> Costs;

	//for (int j = 0; j < data.M; j++) {
	//	StationTotals.push_back(0.0);
	//	Costs.push_back(0.0);
	//}

	bool foundImprovement;
	IloIntArray outlets(env, M);
	//Greedy optimisation
	for (int t = 0; t < data.T; t++) {
		if (Budget[t] < EPS) {
			continue;
		}
		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			StationTotals.clear();
			Costs.clear();
			for (int j = 0; j < data.M; j++) {
				int previous = 0;
				for (int k = 1; k < Mj[j]; k++) {
					if (x_tilde[t][j][k] > 1 - EPS) { previous = k; }
				}
				outlets[j] = previous;
				if (previous < Mj[j] - 1) { StationTotals[make_pair(j, previous + 1)] = 0.0; }
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
						/*for (int k = k0; k < Mj[j]; k++) {
							if (x_tilde[t][j][k-1] > 1 - EPS){
							StationTotals[make_pair(j, k)] += weight;
							break;
							}
						}*/
						if (outlets[j] +1 >= k0) {
								StationTotals[make_pair(j, outlets[j]+1)] += weight;
								break;
							}
					}
				}
			}

			//Reset coverage of stations if can't add new outlet
			for (auto p : StationTotals) {
				int j = p.first.first;
				int k = p.first.second;
				if (k >= Mj[j]) { StationTotals[p.first] = 0.0; continue; }
				double cost = 0.0;
				cost += (double)data.params["c"][t][j];
				if (k == 1) { cost += (double)data.params["f"][t][j]; }
				if (cost > Budget[t]) {
					StationTotals[p.first] = 0.0;
				}
				else {
					//StationTotals[p.first] = p.second / cost;
					Costs[p.first] = cost; 
				}

			}


			pair<int, int> jBar_star = argmax(StationTotals);
			double zBar_star = StationTotals[jBar_star];
			if (zBar_star > 0) {
				//Update budget
				Budget[t] -= Costs[jBar_star];
				int j = jBar_star.first;
				//int k = jBar_star.second;

				for (int tBar = t; tBar < data.T; tBar++) {
					int current = outlets[j];
					//Update solution for current (and all future) years
					x_tilde[tBar][j][current] = 0;
					x_tilde[tBar][j][current+1] = 1;
					//Update coverage of triplets
					for (int i = 0; i < data.N; i++) {
						//double weight = (double)data.params["Ni"][t][i] / (double)R[i];
						for (int r = 0; r < data.R[i]; r++) {
							if ((!coverage[tBar][i][r]) && (data.a[tBar][i][r][j][current+1])) {
								coverage[tBar][i][r] = 1;
							}

						}
					}
				}
				foundImprovement = 1;
			}


		} while (foundImprovement);
	}
	return true;
}

void MulticutCallback::AddPostHeuristic(const IloCplex::Callback::Context& context, IloArray<IloArray<IloNumArray>>& x_tilde)
{
	IloEnv env = context.getEnv();
	
	//Initialise budget
	vector<double> Budget = data.params["B"];
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
	double current = 0.0;
	if (context.getId() == IloCplex::Callback::Context::Id::Relaxation) {
		current = context.getRelaxationValue(budget0);
	}
	else { current = context.getCandidateValue(budget0); }
	//cout << "Spending in year 0: " << current << endl;
	Budget[0] -= current;
	budget0.end();

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
		double current = 0.0;
		if (context.getId() == IloCplex::Callback::Context::Id::Relaxation) {
			current = context.getRelaxationValue(budget);
		}
		else { current = context.getCandidateValue(budget); }
		//cout << "Spending in year " << t << ": " << current << endl;
		Budget[t] -= current;
		budget.end();
	}

	//Initialise coverage
	IloArray<IloArray<IloBoolArray>> coverage(env, T);
	for (int t = 0; t < data.T; t++) {
		//Initialise coverage of triplets
		coverage[t] = IloArray<IloBoolArray>(env, data.N);
		for (int i = 0; i < data.N; i++) {
			coverage[t][i] = IloBoolArray(env, data.R[i]);
			//double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
			for (int r = 0; r < data.R[i]; r++) {
				bool val = 0;
				switch (data.P[t][i][r])
				{
				case triplet::Uncoverable:
					val = 1; //Set to 1 to skip trying to cover later
					break;
				case triplet::Precovered:
					val = 1;
					break;
				default:
					vector<pair<int, int>> cover = data.cover[t][i][r];
					for (long unsigned int jBar = 0; jBar < cover.size(); jBar++) {
						int j = cover[jBar].first;
						int k0 = cover[jBar].second;
						//NOTE: This only uses the integer values from x_tilde deliberately
						for (int k = k0; k < Mj[j]; k++) {
							if (x_tilde[t][j][k] > 1 - EPS) {
								val = 1;
								break;
							}
						}
					}
					break;
				}
				coverage[t][i][r] = val;
			}
		}
	}


	bool feasible = true;
	if (context.getId() == IloCplex::Callback::Context::Id::Relaxation) {
		totalCounter += 1;
		if (heuristic == useHeuristic::WarmstartAndPostSimple) {
			feasible = SimpleRepair(env, x_tilde, Budget, coverage);
			cout << "SimpleRepair complete" << endl;
		}
		else if (heuristic == useHeuristic::WarmstartAndPostGreedy) {
			feasible = GreedyRepair(env, x_tilde, Budget, coverage);
			cout << "GreedyRepair complete" << endl;
		}
		if (feasible) { repairCounter += 1; }

	}
	
	if (feasible) {
		feasible = GreedyFill(env, x_tilde, Budget, coverage);
		cout << "GreedyFill complete" << endl;
	}
	

	if (feasible) {
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
		double obj = 0.0;

		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				//startVar.add(y[t][j]);
				//if (Sol[t][j] > 0) { startVal.add(1); }
				//else { startVal.add(0); }
				int val = 0;
				for (int k = 0; k < Mj[j]; k++) {
					startVar.add(x[t][j][k]);
					startVal.add(x_tilde[t][j][k]);
					if (k >= 1 && x_tilde[t][j][k] > 1 - EPS) { val = 1; }
				}
				startVar.add(y[t][j]);
				startVal.add(val);
			}
		}

		switch (cut_type)
		{
		case multicuts::SingleB1:
		case multicuts::SingleB2:
		{
			double val = 0.0;
			for (int t = 0; t < T; t++) {
					
				for (int i = 0; i < data.N; i++) {
					double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
					for (int r = 0; r < data.R[i]; r++) {
						if ((data.P[t][i][r] == triplet::Precovered) || (data.P[t][i][r] != triplet::Uncoverable && coverage[t][i][r])) { obj += weight; val += weight; }
					}
				}

			}
			startVar.add(theta[0][0]);
			startVal.add(val);
			break;
		}
		case multicuts::Multi1B1:
		case multicuts::Multi1B2:
		{
			for (int t = 0; t < T; t++) {
				double val = 0.0;
				for (int i = 0; i < data.N; i++) {
					double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
					for (int r = 0; r < data.R[i]; r++) {
						if ((data.P[t][i][r] == triplet::Precovered) || (data.P[t][i][r] != triplet::Uncoverable && coverage[t][i][r]) ){ obj += weight; val += weight; }
					}
				}
				startVar.add(theta[t][0]);
				startVal.add(val);
			}
			break;
		}
		case multicuts::Multi2B1:
		case multicuts::Multi2B2:
		{
			for (int i = 0; i < N; i++) {
				double val = 0.0;
				for (int t = 0; t < T; t++) {
					double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
					for (int r = 0; r < data.R[i]; r++) {
						if ((data.P[t][i][r] == triplet::Precovered) || (data.P[t][i][r] != triplet::Uncoverable && coverage[t][i][r])) { obj += weight; val += weight; }
					}
				}
				startVar.add(theta[0][i]);
				startVal.add(val);
			}
			break;
		}
		case multicuts::Multi3B1:
		case multicuts::Multi3B2:
		{
			for (int t = 0; t < T; t++) {
				double val = 0.0;
				for (int i = 0; i < data.N; i++) {
					double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
					for (int r = 0; r < data.R[i]; r++) {
						if ((data.P[t][i][r] == triplet::Precovered) || (data.P[t][i][r] != triplet::Uncoverable && coverage[t][i][r])) { obj += weight; val += weight; }
					}
					startVar.add(theta[t][i]);
					startVal.add(val);
				}
			}
			break;
		}
		default:
			break;
		}


		context.postHeuristicSolution(startVar, startVal, obj, IloCplex::Callback::Context::SolutionStrategy::Propagate);
		cout << "Postheuristic objective: " << obj << endl;
		//for (int j = 0; j < M; j++) {
		//	IloNumArray temp(env, T);
		//	bool post = false;
		//	for (int t = 0; t < T; t++) {
		//		double sum = 0.0;
		//		for (int k = 1; k < Mj[j]; k++) {
		//			sum += k*x_tilde[t][j][k];
		//		}
		//		if (sum > EPS) { temp[t] = sum; post = true; }
		//	}
		//	if (post) { cout << "Station " << j << ": " << temp << endl; }
		//}
	}
	else { cout << "Postheuristic infeasible" << endl; }
	coverage.end();

}

void MulticutCallback::AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde)
{
	IloEnv env = context.getEnv();

	switch (cut_type)
	{
	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B1 cut
	case multicuts::SingleB1:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -EPS) {
				context.rejectCandidate(lhs >= 0).end();
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -EPS) {
				//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
				context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			}
			lhs.end();
			break;
		}
		break; //End of switch statement
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B2 cut
	case multicuts::SingleB2:
	{
		IloExpr lhs(env);
		lhs -= theta[0][0];
		IloNum covered = 0;

		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
		}
		lhs += covered;
		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -EPS) {
				context.rejectCandidate(lhs >= 0).end();
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -EPS) {
				//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
				context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
			}
			lhs.end();
			break;
		}
		break; //End of switch statement
	}
//////////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by year), B1
	case multicuts::Multi1B1:{
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0).end();
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
				}
				lhs.end();
				break;
			}
		}

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by year), B2
	case multicuts::Multi1B2: {
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0).end();
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
				}
				lhs.end();
				break;
			}
		}

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by user class), B1
	case multicuts::Multi2B1: {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < T; t++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0).end();
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
				}
				lhs.end();
				break;
			}
		}

		break; //End of switch statement
	}
	////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by user class), B2
	case multicuts::Multi2B2: {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[0][i];
			IloNum covered = 0;
			for (int t = 0; t < T; t++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0).end();
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
				}
				lhs.end();
				break;
			}
		}

		break; //End of switch statement
	}
//////////////////////////////////////////////////////////////////////////////////////////
//Multi-cut (by year and user class), B1
	case multicuts::Multi3B1: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				switch (context.getId()) {
				case (IloCplex::Callback::Context::Id::Candidate):
					if (context.getCandidateValue(lhs) < -EPS) {
						context.rejectCandidate(lhs >= 0).end();
					}
					lhs.end();
					break;
				case (IloCplex::Callback::Context::Id::Relaxation):
					if (context.getRelaxationValue(lhs) < -EPS) {
						//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
						context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
					}
					lhs.end();
					break;
				}
			}
		}

		break; //End of switch statement
	}
////////////////////////////////////////////////////////////////////////////////////
//Multi-cut (by year and user class), B2
	case multicuts::Multi3B2: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				switch (context.getId()) {
				case (IloCplex::Callback::Context::Id::Candidate):
					if (context.getCandidateValue(lhs) < -EPS) {
						context.rejectCandidate(lhs >= 0).end();
					}
					lhs.end();
					break;
				case (IloCplex::Callback::Context::Id::Relaxation):
					if (context.getRelaxationValue(lhs) < -EPS) {
						//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
						context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse).end();
					}
					lhs.end();
					break;
				}
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
/////////////////////////////////////////////////////////////////////////
}




void MulticutCallback::AddStrengthenedCuts_User(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde, const Num3D& I_tilde){
IloEnv env = context.getEnv();

switch (cut_type)
{
//////////////////////////////////////////////////////////////////////////////////////////
	//Multi-cut (by year), B1
case multicuts::Multi1SB1: {
	for (int t = 0; t < T; t++) {
		IloExpr lhs(env);
		lhs -= theta[t][0];
		IloNum covered = 0;
		for (int i = 0; i < N; i++) {
			IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
			for (int r = 0; r < R[i]; r++) {
				if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
				else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
				else if (data.P[t][i][r] == triplet::Single) {
					if (I_tilde[t][i][r] > 1) {
						covered += weight;
					}
					else {
						int j = data.cover[t][i][r][0].first;
						int k0 = data.cover[t][i][r][0].second;
						//for (int k = k0; k < Mj[j]; k++) {
						//	lhs += weight * x[t][j][k];

						//}
						for (int k = k0; k < Mj[j]-1; k++) {
							if (x_tilde[t][j][k] > EPS){ 
								lhs += weight * x[t][j][k];
								if ((k < Mj[j] - 2) && (x_tilde[t][j][k + 1] < EPS)){ lhs += weight * x[t][j][k+1]; }			
							}
						}
						lhs += weight * x[t][j][Mj[j]-1];
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
							//for (int k = k0; k < Mj[j]; k++) {
							//	lhs += weight * x[t][j][k];

							//}
							for (int k = k0; k < Mj[j] - 1; k++) {
								if (x_tilde[t][j][k] > EPS) {
									lhs += weight * x[t][j][k];
									if ((k < Mj[j] - 2) && (x_tilde[t][j][k + 1] < EPS)) { lhs += weight * x[t][j][k + 1]; }
								}
							}
							lhs += weight * x[t][j][Mj[j] - 1];
						}
					}
				}
			}
		}
		lhs += covered;
		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
			if (context.getCandidateValue(lhs) < -EPS) {
				context.rejectCandidate(lhs >= 0);
			}
			lhs.end();
			break;
		case (IloCplex::Callback::Context::Id::Relaxation):
			if (context.getRelaxationValue(lhs) < -EPS) {
				//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
				context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
			}
			lhs.end();
			break;
		}
	}

	break; //End of switch statement
}
						////////////////////////////////////////////////////////////////////////////////////
						//Multi-cut (by year), B2
case multicuts::Multi1SB2: {
	for (int t = 0; t < T; t++) {
		IloExpr lhs(env);
		lhs -= theta[t][0];
		IloNum covered = 0;
		for (int i = 0; i < N; i++) {
			IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
			for (int r = 0; r < R[i]; r++) {
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
							for (int k = k0; k < Mj[j]; k++) {
								lhs += weight * x[t][j][k];
							}
							//bool coverrest = false;
							//for (int k = k0; k < Mj[j]-1; k++) {
							//	if ((x_tilde[t][j][k] > EPS) || (coverrest)) {
							//		lhs += weight * x[t][j][k];
							//		//coverrest = true;
							//		if ((k < Mj[j] - 2) && (x_tilde[t][j][k + 1] < EPS)) { lhs += weight * x[t][j][k + 1]; }
							//	}
							//}
							//lhs += weight * x[t][j][Mj[j]-1];

						}
					}
				}
			}
		}
		lhs += covered;
		if (context.getRelaxationValue(lhs) < -EPS) {
			//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
			context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse); 
		}
		lhs.end();
	}

	break; //End of switch statement
}
						//////////////////////////////////////////////////////////////////////////////////////////
						//Multi-cut (by year and user class), B1
case multicuts::Multi3B1: {
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[t][i];
			IloNum covered = 0;
			IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
			for (int r = 0; r < R[i]; r++) {
				if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
				else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
				else if (data.P[t][i][r] == triplet::Single) {
					if (I_tilde[t][i][r] > 1) {
						covered += weight;
					}
					else {
						int j = data.cover[t][i][r][0].first;
						int k0 = data.cover[t][i][r][0].second;
						for (int k = k0; k < Mj[j]; k++) {
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
							for (int k = k0; k < Mj[j]; k++) {
								lhs += weight * x[t][j][k];
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0);
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
				}
				lhs.end();
				break;
			}
		}
	}

	break; //End of switch statement
}
						////////////////////////////////////////////////////////////////////////////////////
						//Multi-cut (by year and user class), B2
case multicuts::Multi3B2: {
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < N; i++) {
			IloExpr lhs(env);
			lhs -= theta[t][i];
			IloNum covered = 0;
			IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
			for (int r = 0; r < R[i]; r++) {
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
							for (int k = k0; k < Mj[j]; k++) {
								lhs += weight * x[t][j][k];
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0);
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
				}
				lhs.end();
				break;
			}
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
/////////////////////////////////////////////////////////////////////////
}

void MulticutCallback::AddStrengthenedCuts_Lazy(const IloCplex::Callback::Context& context, const IloArray<IloNumArray>& x_tilde, const Num3D& I_tilde)
{
	IloEnv env = context.getEnv();

	switch (cut_type)
	{
		//////////////////////////////////////////////////////////////////////////////////////////
			//Multi-cut (by year), B1
	case multicuts::Multi1SB1: {
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							//for (int k = k0; k < Mj[j]; k++) {
							//	lhs += weight * x[t][j][k];
							//}
							int kBar = max((int) x_tilde[t][j], k0);
							for (int k = kBar; k<Mj[j];k++) {
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
								//for (int k = k0; k < Mj[j]; k++) {
								//	lhs += weight * x[t][j][k];
								//}
								int kBar = max((int)x_tilde[t][j], k0);
								for (int k = kBar; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
			}
			lhs += covered;
			switch (context.getId()) {
			case (IloCplex::Callback::Context::Id::Candidate):
				if (context.getCandidateValue(lhs) < -EPS) {
					context.rejectCandidate(lhs >= 0);
				}
				lhs.end();
				break;
			case (IloCplex::Callback::Context::Id::Relaxation):
				if (context.getRelaxationValue(lhs) < -EPS) {
					//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
					context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
				}
				lhs.end();
				break;
			}
		}

		break; //End of switch statement
	}
							 ////////////////////////////////////////////////////////////////////////////////////
							 //Multi-cut (by year), B2
	case multicuts::Multi1SB2: {
		for (int t = 0; t < T; t++) {
			IloExpr lhs(env);
			lhs -= theta[t][0];
			IloNum covered = 0;
			for (int i = 0; i < N; i++) {
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
								//int kBar = max((int)x_tilde[t][j], k0);
								//for (int k = kBar; k < Mj[j]; k++) {
								//	lhs += weight * x[t][j][k];
								//}
							}
						}
					}
				}
			}
			lhs += covered;
			if (context.getCandidateValue(lhs) < -EPS) {
				context.rejectCandidate(lhs >= 0);
			}
			lhs.end();
		}

		break; //End of switch statement
	}
							//////////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year and user class), B1
	case multicuts::Multi3B1: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (data.P[t][i][r] == triplet::Uncoverable) { continue; }
					else if (data.P[t][i][r] == triplet::Precovered) { covered += weight; }
					else if (data.P[t][i][r] == triplet::Single) {
						if (I_tilde[t][i][r] > 1) {
							covered += weight;
						}
						else {
							int j = data.cover[t][i][r][0].first;
							int k0 = data.cover[t][i][r][0].second;
							for (int k = k0; k < Mj[j]; k++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				switch (context.getId()) {
				case (IloCplex::Callback::Context::Id::Candidate):
					if (context.getCandidateValue(lhs) < -EPS) {
						context.rejectCandidate(lhs >= 0);
					}
					lhs.end();
					break;
				case (IloCplex::Callback::Context::Id::Relaxation):
					if (context.getRelaxationValue(lhs) < -EPS) {
						//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
						context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
					}
					lhs.end();
					break;
				}
			}
		}

		break; //End of switch statement
	}
							////////////////////////////////////////////////////////////////////////////////////
							//Multi-cut (by year and user class), B2
	case multicuts::Multi3B2: {
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < N; i++) {
				IloExpr lhs(env);
				lhs -= theta[t][i];
				IloNum covered = 0;
				IloNum weight = (IloNum)data.params["Ni"][t][i] / (IloNum)R[i];
				for (int r = 0; r < R[i]; r++) {
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
								for (int k = k0; k < Mj[j]; k++) {
									lhs += weight * x[t][j][k];
								}
							}
						}
					}
				}
				lhs += covered;
				switch (context.getId()) {
				case (IloCplex::Callback::Context::Id::Candidate):
					if (context.getCandidateValue(lhs) < -EPS) {
						context.rejectCandidate(lhs >= 0);
					}
					lhs.end();
					break;
				case (IloCplex::Callback::Context::Id::Relaxation):
					if (context.getRelaxationValue(lhs) < -EPS) {
						//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse);
						context.addUserCut(lhs >= 0, IloCplex::UseCutPurge, IloFalse);
					}
					lhs.end();
					break;
				}
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
	/////////////////////////////////////////////////////////////////////////
}

