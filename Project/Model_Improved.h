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
	double ObjectiveValue = -1;
	double SolveTime = -1;
	double TotalTime = -1;
	ArrayXXd Solution;
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


