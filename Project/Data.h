#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
#include <exception>
#include <map>

#include "Librairies/nlohmann/json.hpp"
#include "Librairies/Eigen/Core"
#include "Librairies/Eigen/SparseCore"

#define EPS 1.0e-6

using json = nlohmann::json;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using namespace std;

typedef Eigen::SparseMatrix<double> SparseXXd;
typedef Eigen::SparseMatrix<bool> SparseXXb;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;
typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::Triplet<bool> Tripletb;


class Data
{

public:
	//Number of time periods
	int T;
	//Number of station-outlet pairs, time period invariant
	int M_bar;
	//Number of (reduced) triplets, divided by time period
	vector<int> P;
	//Parameters from the Shared.json file
	json params;

	//In the data, all stations are assumed to be level 3, but this parameter can be adjusted if desired
	//(only usable when loading from from utilities) 
	vector<bool> isLevel3 = { 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0 };

	//Precomputed coverage of each station-outlet pair, for each t
	vector<SparseXXd> a;

	//Precomputed values for a_{ij}^t * d_{ij}^t, for each t
	vector<SparseXXd> CutCoeffs;

	//Total demand covered *only* by each station-outlet pair, for each t
	//Called "Js" in the paper
	vector<Eigen::VectorXd> Ps;

	//Total demand that is covered by home charging, for each t
	vector<double> Precovered;

	//Vector of weights d_{ij}^t, for each t
	vector<Eigen::VectorXd> weights;

	/*
	Load parameters from 'sharedFile' and coverage for 'instanceFile'. 
	This method is more efficient than load_fromUtilities as the coverage is already computed
	and must simply be loaded in. In general, this method should be used if you are not 
	changing the data.
	*/
	void load(string sharedFile, string instanceFile, bool verbose);

	/*
	Load parameters from 'sharedFile' and compute coverage for 'instanceFile'.
	This method requires calculating the utilities and resulting coverage values, 
	and so should only be used if you need to change the underlying utilities
	(e.g. if using level 2 and level 3 chargers and adjusting costs/utilities)
	*/
	void load_fromUtilities(string sharedFile, string instanceFile, bool verbose = false);
	
	/*
	Deprecated function for calculating the total demand covered by a solution. Solutions
	are represented as Arrays instead of as vectors. 
	*/
	double SolutionQuality(vector<vector<int>> Sol);

	/*
	Function for calculating the total demand covered by a solution.
	*/
	double SolutionQuality(ArrayXXd Sol);

private:
	/*
	Loads the coverage for an instance file. Only used as part of load function, and
	helps clean up the code a bit.
	*/
	void create_covering(string instanceFile, bool verbose);


};

