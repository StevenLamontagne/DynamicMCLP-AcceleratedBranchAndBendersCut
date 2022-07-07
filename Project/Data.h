#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>

#ifdef _WIN32
#include <nlohmann/json.hpp>
#endif

#ifdef __linux__
#include "json.hpp"
#endif


using json = nlohmann::json;
using namespace std;

enum class triplet:char {
	Uncoverable,
	Precovered,
	Single,
	Multi
};

class Data
{
public:
	int T;
	int N;
	int M;
	vector<int> Mj;
	vector<int> R;

	json params;
	vector<vector<vector<vector<vector<bool>>>>> a;
	vector<vector<vector<bool>>> home;
	vector<vector<vector<triplet>>> P;
	vector<vector<vector<vector<pair<int, int>>>>> cover;

	void load(string file, bool verbose);
	double SolutionQuality(vector<vector<int>> Sol);

private:
	void create_covering(bool verbose);
};

