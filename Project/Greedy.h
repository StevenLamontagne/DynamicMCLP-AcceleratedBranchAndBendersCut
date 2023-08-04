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


/*
Class for storing and solving the dynamic MCLP (EV charging station placement problem) using the
Greedy method from Lamontagne et al. (2022). This is used as a warmstart for A-B&BC and A-B&BC+LB methods
but can also be used as a standalone approach.

'budgetType' is an option for using either knapsack-style or cardinality-style budget constraints. This was not
in the final paper, but was left in case it could be useful reference.
*/
class Greedy
{
public:
	void SetData(const Data& newData) { data = newData; };
	void Solve(bool _verbose = false, BUDGET_TYPE budgetType = BUDGET_TYPE::Knapsack);

	//Storing solving statistics
	double SolveTime;
	vector<vector<int>> Solution;
	double SolutionQuality = 0.0;

	bool verbose = false;

private:
	Data data;

	//These are the internal equivalents for the coverage and budget found in Data class
	//since the greedy methods updates the values as outlets are added
	vector<vector<vector<bool>>> coverage;	
	vector<double> Budget;

	/*
	Simple function for returning the index of the maximum element. Used for
	identifying best outlet to place.
	*/
	int argmax(vector<double> vec);



};

