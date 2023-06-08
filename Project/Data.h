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

enum class priceProfile
{
	MixedLevel2andLevel3_Unperturbed,
	MixedLevel2andLevel3_Perturbed,
	OnlyLevel3_Unperturbed,
	OnlyLevel3_Perturbed
};

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

	vector<bool> isLevel3 = { 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0 };
	vector<float> costMultiplier = { 0.7707, 0.9917, 1.1381, 0.803, 1.1038, 0.9992, 1.1902, 1.1833, 1.2401, 1.1188, 0.8373,
									 1.0819, 1.0997, 0.8856, 1.0984, 0.8053, 0.8103, 1.172, 0.7956, 1.1326, 1.1701, 1.1638,
									 0.7603, 0.8443, 0.9484, 1.0118, 0.9439, 0.9442, 0.9484, 0.8184 };


	vector<SparseXXd> a;
	vector<SparseXXd> CutCoeffs;
	vector<Eigen::VectorXd> Ps;
	vector<double> Precovered;
	vector<Eigen::VectorXd> weights;

	void load(string sharedFile, string instanceFile, bool verbose);
	void load_fromUtilities(string sharedFile, string instanceFile, bool verbose = false, priceProfile priceProfile = priceProfile::MixedLevel2andLevel3_Perturbed);
	
	double SolutionQuality(vector<vector<int>> Sol);
	double SolutionQuality(ArrayXXd Sol);

private:
	void create_covering(string instanceFile, bool verbose);


};

