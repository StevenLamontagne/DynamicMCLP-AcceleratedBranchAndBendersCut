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
	void Solve(json params);

	void GetSolution(IloCplex& cplex, BoolVar3D& x);

	//Storing information about solution
	float ObjectiveValue = -1;
	float SolveTime = -1;
	float TotalTime = -1;
	vector<vector<int>> Solution;
	double OptimalityGap = -1;

	//Options
	multicuts multicut = multicuts::Multi1B1;
	useHeuristic heuristic = useHeuristic::Warmstart;
	int nGRASP = 0;
	bool verbose = false;
	double overlap_threshold = 0.1;
	double grasp_threshold = 0.85;
	bool use_trust = false;
	double trust_threshold = 2.5;

	//Storing solving statistics
	map<string, int> stats;



	

private:
	Data data;

	//vector<map<pair<int, int>, vector<pair<int, int>>>> cover_station;
	vector<vector<vector<vector<pair<int, int>>>>> cover_station;
	vector<map<pair<int, int>, map<pair<int, int>, double>>> overlap;

	bool GRASPCut(IloEnv& env, IloModel& model, BoolVar3D& x, IloArray<IloNumVarArray>& theta);

	


};

