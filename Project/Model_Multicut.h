#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

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

private:
	Data data;


};

