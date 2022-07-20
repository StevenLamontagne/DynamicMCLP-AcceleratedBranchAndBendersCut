#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

#include "Data.h"
#include "CoverageCallback.h"

class Model_BendersCordeau
{
public:
	void SetData(const Data & newData) { data = newData; };
	void Solve(cuts cut_type, vector<vector<int>>* warmstart = nullptr);

	float ObjectiveValue;
	float SolveTime;
	vector<vector<int>> Solution;
	float OptimalityGap;


private:
	Data data;


};

