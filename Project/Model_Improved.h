#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <cmath>

#include "Data_Improved.h"
#include "Callback_Improved.h"
#include "Greedy_Improved.h"

class Model_Improved
{
public:
	void SetData(const Data_Improved& newData);
	void Solve(json params);

	void GetSolution(IloCplex& cplex, BoolVar2D& x);

	//Storing information about solution
	float ObjectiveValue = -1;
	float SolveTime = -1;
	float TotalTime = -1;
	vector<vector<int>> Solution;
	double OptimalityGap = -1;

	//Options
	bool verbose = false;
	double overlap_threshold = 0.1;
	bool use_trust = false;
	double trust_threshold = 2.0;

	//Storing solving statistics
	map<string, int> stats;





private:
	Data_Improved data;





};


