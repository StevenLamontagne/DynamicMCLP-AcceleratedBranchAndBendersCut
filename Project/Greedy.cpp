#include "Greedy.h"

void Greedy::Solve(bool _verbose, BUDGET_TYPE budgetType)
{
	verbose = _verbose;

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	int M = data.params["M"];
	vector<int> Mj = data.params["Mj"];

	map<pair<int,int>, int> coordConverter;
	for (int j_bar = 0; j_bar < M_bar; j_bar++) { coordConverter[make_pair(data.params["station_coord"][j_bar][0], data.params["station_coord"][j_bar][1])] = j_bar; }

	vector<VectorXd> cover;
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
	//Initialisation
	for (int t = 0; t < T; t++) {
		//Initialise solution to x0
		vector<int>  Sol1;
		for (int j = 0; j < M; j++) {
			int val = 0;
			for (int k = 1; k < Mj[j]; k++) {
				if (data.params["x0"][j][k] == 1) { val = k; }
			}
			Sol1.push_back(val);
		}
		Solution.push_back(Sol1);

		//Initialise budget
		Budget.push_back(data.params["B"][t]);

		//Initialise coverage of triplets and solution quality
		cover.push_back(VectorXd::Constant(P[t], 0.0)); //By construction, all of the precovered triplets have been removed
		//SolutionQuality += data.Precovered[t];
		
	}
	if (verbose) {
		cout << "Initial coverage: " << SolutionQuality << endl << "\n" << endl;
	}


	//Initialise vectors for holding information for each station
	VectorXd StationTotals = VectorXd::Constant(M, 0.0);
	VectorXd Costs = VectorXd::Constant(M, 0.0);

	//Greedy optimisation
	for (int t = 0; t < T; t++) {
		if (verbose) { cout << "Starting year " << t << endl; }

		bool foundImprovement = 0;
		
		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			StationTotals.setConstant(0.0);
			Costs.setConstant(0.0);
			VectorXd uncovered = (cover[t].array() < 1).matrix().cast<double>();
			//cout << uncovered.sum() << endl;;
			
			int StationsOpened = 0;
			int OutletsInstalled = 0;
			for (int j = 0; j < M; j++) {
				if (Solution[t][j] == 0) { continue; }
				if (t == 0) {
					if (Solution[t][j] > 0) { StationsOpened += 1; }
					OutletsInstalled += Solution[t][j];
				}
				else {
					if (Solution[t - 1][j] == 0) { StationsOpened += 1; }
					OutletsInstalled += Solution[t][j] - Solution[t-1][j];
				}
				
			}

			//Iterate over stations
			for (int j = 0; j < M; j++) {
				int updated = Solution[t][j] + 1;
				//Skip: Maximum number of outlets already placed
				if (updated >= Mj[j]) {
					//if (verbose) { cout << "Skipping station " << j << " due to maximum outlets." << endl; }
					continue;
				}
				double cost = data.params["c"][t][j][updated];

				//Skip:Budget insufficient for new outlet
				if ((budgetType == BUDGET_TYPE::Knapsack) && (cost > Budget[t])){	
					continue;
				}	
				else if (budgetType == BUDGET_TYPE::OutletCount) {
					//Skip:Too many outlets installed or stations opened
					if (t == 0) {
						if ((Solution[t][j] == 0) && (StationsOpened >= data.params["Stations_maxNewStationsPerTimePeriod"])) { continue; }
						if (OutletsInstalled >= data.params["Stations_maxNewOutletsPerTimePeriod"]) { continue; }
					}
					else {
						if ((Solution[t - 1][j] == 0) && (Solution[t][j] == 0) && (StationsOpened >= data.params["Stations_maxNewStationsPerTimePeriod"])) { continue; }
						if (OutletsInstalled >= data.params["Stations_maxNewOutletsPerTimePeriod"]) { continue; }
					}
				}
				Costs(j) = cost;

				int j_bar = coordConverter[make_pair(j, updated)];
				StationTotals(j) = data.CutCoeffs[t].col(j_bar).dot(uncovered);
				StationTotals(j) += data.Ps[t](j_bar);
				if (budgetType == BUDGET_TYPE::Knapsack) { StationTotals(j) /= cost; }  //Should we care about the cost in the number of outlets case? 			
			}

			//cout << "Station totals: [";
			//for (double val : StationTotals) { cout << val << ", "; }
			//cout << "]" << endl;

			int j_star;
			double z_star = StationTotals.maxCoeff(&j_star);
			if (z_star > 0) {
				if (verbose) {
					cout << "Found improvement of amount " << z_star << endl;
					cout << "Adding outlet at station " << j_star << endl;
				}
				//Update budget
				Budget[t] -= Costs(j_star);

				for (int tBar = t; tBar < T; tBar++) {
					//Update solution for current (and all future) years
					int updated = Solution[tBar][j_star] + 1;
					Solution[tBar][j_star] += 1;
					cover[tBar] += data.a[tBar].col(coordConverter[make_pair(j_star, updated)]);
				}
				foundImprovement = 1;
			}


		} while (foundImprovement);
		if (verbose) { cout << "\n" << endl; }
	}
	SolveTime = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
	SolutionQuality = data.SolutionQuality(Solution);
}


int Greedy::argmax(vector<double> vec)
{
	auto maxVal = max_element(vec.begin(), vec.end());
	int argmaxVal = distance(vec.begin(), maxVal);
	return argmaxVal;
}
