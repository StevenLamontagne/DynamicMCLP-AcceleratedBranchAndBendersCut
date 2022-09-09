#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
#include <exception>

#ifdef _WIN32
#include <nlohmann/json.hpp>
#endif

#ifdef __linux__
//# define JSON_DIAGNOSTICS 1
#include "json.hpp"
#endif


using json = nlohmann::json;
using namespace std;

enum class triplet :char {
	Uncoverable,
	Precovered,
	Single,
	Multi
};

class Data_Strengthened
{
public:
	int T;
	int N;
	int M;
	vector<int> Mj;
	vector<int> R;

	json params;
	vector<vector<vector<vector<vector<bool>>>>> aBar;
	vector<vector<vector<bool>>> home;
	vector<vector<vector<triplet>>> P;
	vector<vector<vector<vector<pair<int, int>>>>> cover;

	void load(string file, bool verbose);
	double SolutionQuality(vector<vector<int>> Sol);

	template <class Temp>
	vector<double> RemainingBudget(const vector<Temp>& Sol) {
		vector<double> newBudget = params["B"];
		for (int t = 0; t < T; t++) {
			for (int j = 0; j < M; j++) {
				double previous = 0.0;
				if (t == 0) { previous = params["x0"][j]; }
				else { previous = Sol[t - 1][j]; }

				double delta = Sol[t][j] - previous;
				if (previous < 1) { newBudget[t] -= (double)params["f"][t][j] * delta; }
				newBudget[t] -= (double)params["c"][t][j] * delta;
			}
		}
		return newBudget;
	};

private:
	void create_covering(bool verbose);
};
