#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

#include "Data.h"

class Model_BendersCordeau
{
public:
	void SetData(const Data & newData) { data = newData; };
	void Solve(vector<vector<int>>* warmstart = nullptr);

	float ObjectiveValue;
	float SolveTime;
	vector<vector<int>> Solution;


private:
	Data data;


};

