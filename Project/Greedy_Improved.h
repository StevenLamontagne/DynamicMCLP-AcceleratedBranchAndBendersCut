#pragma once
#include "Data_Improved.h"

#include <exception>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>

using Eigen::VectorXd;
typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayVectord;

class Greedy_Improved
{
public:
	void SetData(const Data_Improved& newData) { data = newData; };
	void Solve(bool _verbose = false);

	double SolveTime;
	vector<vector<int>> Solution;
	double SolutionQuality = 0.0;

	bool verbose = false;

private:
	Data_Improved data;
	vector<vector<vector<bool>>> coverage;
	int argmax(vector<double> vec);
	vector<double> Budget;



};

