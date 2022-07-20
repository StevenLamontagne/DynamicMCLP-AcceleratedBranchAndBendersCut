#include "Greedy.h"

void Greedy::Solve()
{
	//Define constants for easier reading and writing
	int& T = data.T;
	int& M = data.M;
	int& N = data.N;
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;

	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
	//Initialisation
	for (int t = 0; t < T; t++) {
		//Initialise solution to x0
		vector<int>  Sol1;
		for (int j = 0; j < M; j++) {
			Sol1.push_back(data.params["x0"][j]);
		}
		Solution.push_back(Sol1);

		//Initialise budget
		Budget.push_back(data.params["B"][t]);

		//Initialise coverage of triplets
		vector<vector<bool>> cover1;
		for (int i = 0; i < N; i++) {
			vector<bool> cover2;
			double weight = (double)data.params["Ni"][t][i] / (double)R[i];
			for (int r = 0; r < R[i]; r++) {
				bool val = 0;
				switch (data.P[t][i][r])
				{
				case triplet::Uncoverable:
					val = 1; //Set to 1 to skip trying to cover later
					break;
				case triplet::Precovered:
					val = 1;
					SolutionQuality += weight;
					break;
				default:
					break;
				}
				cover2.push_back(val);
			}
			cover1.push_back(cover2);
		}
		coverage.push_back(cover1);
	}
	cout << "Initial coverage: " << SolutionQuality << endl;
	cout << "\n" << endl;

	//Initialise vectors for holding information for each station
	vector<double> StationTotals;
	vector<double> Costs;
	for (int j = 0; j < M; j++) {
		StationTotals.push_back(0.0);
		Costs.push_back(0.0);
	}

	//Greedy optimisation
	for (int t = 0; t < T; t++) {
		cout << "Starting year " << t << endl;
		bool foundImprovement;
		
		//Check for coverage from initial solution/solution from last year
		for (int i = 0; i < N; i++) {
			double weight = (double)data.params["Ni"][t][i] / (double)R[i];
			for (int r = 0; r < R[i]; r++) {
				if (coverage[t][i][r]) { continue; }
				vector<pair<int, int>> cover = data.cover[t][i][r];
				for (int jBar = 0; jBar < cover.size(); jBar++) {
					int j = cover[jBar].first;
					int k0 = cover[jBar].second;
					if (Solution[t][j] >= k0) {
						coverage[t][i][r] = 1;
					}
				}
				if (coverage[t][i][r]) { SolutionQuality += weight; }
			}
		}

		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			double coverable = 0.0;
			for (int j = 0; j < M; j++) {
				StationTotals[j] = 0.0;
				Costs[j] = 0.0;
			}

			//Check for coverage from new stations
			for (int i = 0; i < N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (coverage[t][i][r]) { continue; }
					vector<pair<int, int>> cover = data.cover[t][i][r];
					for (int jBar = 0; jBar < cover.size(); jBar++) {
						int j = cover[jBar].first;
						int k0 = cover[jBar].second;
						if (Solution[t][j] + 1 >= k0) {
								StationTotals[j] += weight;
						}
					}
				}
			}
	
			//Reset coverage of stations to 0 if can't add new outlet
			for (int j = 0; j < M; j++) {
				//Maximum number of outlets already placed
				if (Solution[t][j] >= Mj[j] - 1) { StationTotals[j] = 0.0; }

				//Budget insufficient for new outlet
				Costs[j] += data.params["c"][t][j];;
				if (Solution[t][j] == 0) { Costs[j] += data.params["f"][t][j]; }
				if (Costs[j] > Budget[t]) {
					StationTotals[j] = 0.0;
				}
				//cout << "Coverage of station " << j << ": " << StationTotals[j] << endl;
			}


			int j_star = argmax(StationTotals);
			double z_star = StationTotals[j_star];
			if (z_star > 0) {
				cout << "Found improvement of amount " << z_star << endl;
				cout << "Adding outlet at station " << j_star << endl;
				//Update solution for current (and all future) years
				for (int tBar = t; tBar < T; tBar++) {
					Solution[tBar][j_star] += 1;
				}

				//Update budget
				Budget[t] -= Costs[j_star];

				//Update solution quality (Solution quality unnecessary?)
				SolutionQuality += z_star;

				//Update coverage of triplets
				for (int i = 0; i < N; i++) {
					for (int r = 0; r < R[i]; r++) {
						if (!coverage[t][i][r]) {
							//Check for pair directly
							if (std::find(data.cover[t][i][r].begin(), data.cover[t][i][r].end(), std::make_pair(j_star, Solution[t][j_star])) != data.cover[t][i][r].end()) {
								coverage[t][i][r] = 1;
							}
							//else {
							//	coverable += (double)data.params["Ni"][t][i] / (double)R[i];
							//}
						}

						////Check for coverage more manually
						//vector<pair<int, int>> cover = data.cover[t][i][r];
						//for (int jBar = 0; jBar < cover.size(); jBar++) {
						//	int j = cover[jBar].first;
						//	int k0 = cover[jBar].second;
						//	if (Solution[t][j] >= k0){
						//		coverage[t][i][r] = 1;
						//	}
						//}
					}
				}
				foundImprovement = 1;	
				//cout << "Coverable left: " << coverable << endl;
			}


		} while (foundImprovement);
		cout << "\n" << endl;
	}
	SolveTime = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
}


int Greedy::argmax(vector<double> vec)
{
	auto maxVal = max_element(vec.begin(), vec.end());
	int argmaxVal = distance(vec.begin(), maxVal);
	return argmaxVal;
}
