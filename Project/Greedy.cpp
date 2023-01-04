#include "Greedy.h"

void Greedy::Solve(bool _verbose)
{	
	verbose = _verbose;

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
			int val = 0;
			for (int k = 1; k < Mj[j]; k++) {
				if (data.params["x0"][j] == 1) { val = k; }
			}
			Sol1.push_back(val);
		}
		Solution.push_back(Sol1);

		//Initialise budget
		Budget.push_back(data.params["B"][t]);

		//Initialise coverage of triplets
		vector<vector<bool>> cover1;
		for (int i = 0; i < N; i++) {
			vector<bool> cover2;
			for (int r = 0; r < R[i]; r++) {
				bool val = 0;
				switch (data.P[t][i][r])
				{
				case triplet::Uncoverable:
					val = 1; //Set to 1 to skip trying to cover_triplet later
					break;
				case triplet::Precovered:
					val = 1;
					break;
				default:
					for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
						int j = cover_triplet.first;
						int k0 = cover_triplet.second;
						if (Solution[t][j] >= k0) {
							val = 1;
						}
					}
					break;
				}
				cover2.push_back(val);
			}
			cover1.push_back(cover2);
		}
		coverage.push_back(cover1);
	}
	if (verbose) {
		cout << "Initial coverage: " << SolutionQuality << endl << "\n" << endl;
	}
	

	//Initialise vectors for holding information for each station
	vector<double> StationTotals;
	vector<double> Costs;
	for (int j = 0; j < M; j++) {
		StationTotals.push_back(0.0);
		Costs.push_back(0.0);
	}

	//Greedy optimisation
	for (int t = 0; t < T; t++) {
		if (verbose) { cout << "Starting year " << t << endl; }
		bool foundImprovement;
		
		//Start adding outlets
		do
		{
			//Reset values
			foundImprovement = 0;
			for (int j = 0; j < M; j++) {
				StationTotals[j] = 0.0;
				Costs[j] = 0.0;
			}

			//Check for coverage by adding outlet
			for (int i = 0; i < N; i++) {
				double weight = (double)data.params["Ni"][t][i] / (double)R[i];
				for (int r = 0; r < R[i]; r++) {
					if (coverage[t][i][r]) { continue; }
					for (pair<int, int> cover_triplet : data.cover_triplet[t][i][r]) {
						int j = cover_triplet.first;
						int k0 = cover_triplet.second;
						if (Solution[t][j] + 1 >= k0) {
							StationTotals[j] += weight;
						}
					}
				}
			}
	
			//Reset coverage of stations to 0 if can't add new outlet
			for (int j = 0; j < M; j++) {
				//Maximum number of outlets already placed
				if (Solution[t][j] >= Mj[j] - 1) { 
					StationTotals[j] = 0.0;
					continue;
				}

				//Budget insufficient for new outlet
				Costs[j] = (double) data.params["c"][t][j][Solution[t][j]+1];
				if (Costs[j] > Budget[t]) {
					StationTotals[j] = 0.0;
				}
			}


			int j_star = argmax(StationTotals);
			double z_star = StationTotals[j_star];
			if (z_star > 0) {
				if (verbose) {
					cout << "Found improvement of amount " << z_star << endl;
					cout << "Adding outlet at station " << j_star << endl;
				}
				//Update budget
				Budget[t] -= Costs[j_star];

				for (int tBar = t; tBar < T; tBar++) {
					//Update solution for current (and all future) years
					Solution[tBar][j_star] += 1;

					//Update coverage of triplets
					for (int i = 0; i < N; i++) {
						for (int r = 0; r < R[i]; r++) {
							if ((!coverage[tBar][i][r]) && (data.a[tBar][i][r][j_star][Solution[tBar][j_star]])) {
								coverage[tBar][i][r] = 1;								
							}
						}
					}
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
