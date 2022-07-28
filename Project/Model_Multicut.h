#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <cmath>

#include "Data.h"
#include "MulticutCallback.h"

class Model_Multicut
{
public:
	void SetData(const Data& newData) { data = newData; };
	void Solve(multicuts cut_type, vector<vector<int>>* warmstart = nullptr);

	float ObjectiveValue;
	float SolveTime;
	vector<vector<int>> Solution;
	float OptimalityGap;

	bool GreedyRepair(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage);
	bool GreedyFill(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage);


private:
	Data data;

	int argmax(vector<double> vec) {
		auto maxVal = max_element(vec.begin(), vec.end());
		int argmaxVal = distance(vec.begin(), maxVal);
		return argmaxVal;
	};

	template <class T>
	double sum(T vec) {
		double val = 0.0;
		for (auto i : vec) { val += i; }
		return val;
	};

	vector<int> GetFractional(vector<vector<double>> Sol);


};

