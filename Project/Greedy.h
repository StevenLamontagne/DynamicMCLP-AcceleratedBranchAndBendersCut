#pragma once
#include "Data.h"
#include "Utils.h"

#include <exception>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>

using Eigen::VectorXd;
typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayVectord;



class Greedy
{
public:
	void SetData(const Data& newData) { data = newData; };
	void Solve(bool _verbose = false, BUDGET_TYPE budgetType = BUDGET_TYPE::Knapsack);

	double SolveTime;
	vector<vector<int>> Solution;
	double SolutionQuality = 0.0;

	bool verbose = false;

private:
	Data data;
	vector<vector<vector<bool>>> coverage;
	int argmax(vector<double> vec);
	vector<double> Budget;



};

