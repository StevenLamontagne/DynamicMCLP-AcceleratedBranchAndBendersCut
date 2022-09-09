#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <cmath>

#include "Data_Strengthened.h"
#include "MulticutCallback_Strengthened.h"
//#include "Greedy.h"

class Model_Strengthened
{
public:
	void SetData(const Data_Strengthened& newData) { data = newData; };
	void Solve(multicuts cut_type, useHeuristic heuristic, int nGRASP);

	float ObjectiveValue;
	float SolveTime;
	float TotalTime;
	vector<vector<int>> Solution;
	float OptimalityGap;
	multicuts multicut = multicuts::Multi1B1;

	double RepairRatio = 0;

	bool GreedyRepair(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage);
	bool GreedyFill(vector<vector<double>>& Sol, vector<double>& Budget, vector<vector<vector<bool>>>& coverage);
	bool GRASPCut(IloEnv& env, IloModel& model, BoolVar3D& x, IloArray<IloNumVarArray>& theta);

private:
	Data_Strengthened data;

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

