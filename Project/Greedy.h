#pragma once
#include "Data.h"

#include <exception>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>




class Greedy
{
public:
	void SetData(const Data& newData) { data = newData; };
	void Solve();

	double SolveTime;
	vector<vector<int>> Solution;
	double SolutionQuality = 0.0;
	
private:
	Data data;
	vector<vector<vector<bool>>> coverage;
	int argmax(vector<double> vec);
	vector<double> Budget;




};

