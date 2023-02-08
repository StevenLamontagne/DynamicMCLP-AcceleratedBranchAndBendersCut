#pragma once

#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

#include "Data.h"



class Model_BB
{
public:
	void SetData(const Data & newData) { data = newData; };
	void Solve();

	float ObjectiveValue = -1;
	float SolveTime = -1;
	float OptimalityGap = -1;
	int nNodes = 0;
	vector<vector<int>> Solution;


private: 
	Data data;

};

