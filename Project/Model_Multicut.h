#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <cmath>

#include "Data.h"
#include "MulticutCallback.h"
#include "Greedy.h"

class Model_Multicut
{
public:
	void SetData(const Data& newData);
	void Solve(multicuts _multicut = multicuts::Multi1B1, useHeuristic _heuristic = useHeuristic::Warmstart, int _nGRASP = 0, bool _verbose = false);

	//Storing information about solution
	float ObjectiveValue;
	float SolveTime;
	float TotalTime;
	vector<vector<int>> Solution;
	float OptimalityGap;

	//Options
	multicuts multicut = multicuts::Multi1B1;
	useHeuristic heuristic = useHeuristic::Warmstart;
	int nGRASP = 0;
	bool verbose = false;
	double threshold = 0.1;

	//Storing solving statistics
	map<string, int> stats;



	

private:
	Data data;

	//vector<map<pair<int, int>, vector<pair<int, int>>>> cover_station;
	vector<vector<vector<vector<pair<int, int>>>>> cover_station;
	vector<map<pair<int, int>, map<pair<int, int>, double>>> overlap;

	bool GRASPCut(IloEnv& env, IloModel& model, BoolVar3D& x, IloArray<IloNumVarArray>& theta);

	


};

