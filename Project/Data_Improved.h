#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
#include <exception>
#include <map>

#include "Librairies/nlohmann/json.hpp"
#include "Librairies/Eigen/Core"
#include "Librairies/Eigen/SparseCore"

#define EPS 1.0e-6

using json = nlohmann::json;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using namespace std;

typedef Eigen::SparseMatrix<double> SparseXXd;
typedef Eigen::SparseMatrix<bool> SparseXXb;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;
typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::Triplet<bool> Tripletb;



class Data_Improved
{

public:
	//Number of time periods
	int T;
	//Number of station-outlet pairs, time period invariant
	int M_bar;
	//Number of (reduced) triplets, divided by time period
	vector<int> P;
	//Parameters from the Shared.json file
	json params;


	vector<SparseXXd> a;
	vector<SparseXXd> CutCoeffs;
	vector<Eigen::VectorXd> Ps;
	vector<double> Precovered;
	vector<Eigen::VectorXd> weights;

	void load(string set, string file, bool verbose);
	double SolutionQuality(vector<vector<int>> Sol);

	template <class Temp>
	vector<double> RemainingBudget(const vector<Temp>& Sol) {
		vector<double> newBudget = params["B"];
		//for (int t = 0; t < T; t++) {
		//	for (int j = 0; j < M; j++) {
		//		double previous = 0.0;
		//		if (t == 0) { previous = params["x0"][j]; }
		//		else { previous = Sol[t - 1][j]; }

		//		double delta = Sol[t][j] - previous;
		//		if (previous < 1) { newBudget[t] -= (double)params["f"][t][j] * delta; }
		//		newBudget[t] -= (double)params["c"][t][j] * delta;
		//	}
		//}
		return newBudget;
	};

private:
	void create_covering(string path, string test, bool verbose);
};

