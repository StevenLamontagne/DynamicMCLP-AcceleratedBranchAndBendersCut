#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <cmath>
#include <stdexcept>

#include "Data_Improved.h"
#include "Callback_LocalBranching.h"
#include "Greedy_Improved.h"

class Model_LocalBranching
{
public:
	void SetData(const Data_Improved& newData);
	void Solve(json params);

	void GetSolution(IloCplex& cplex, BoolVar2D& x);

	~Model_LocalBranching() {
		model.end();
		env.end();
	}
	

	bool verbose = false;
	//Storing solving statistics
	map<string, int> stats;


private:
	//Define variables and references
	Data_Improved data;
	IloEnv env;
	IloModel model;
	int T;
	int M_bar;
	vector<int> P;
	BoolVar2D x;
	IloNumVarArray theta;
	IloNumVar theta_obj;
	BoolVar2D trust;

	//Storing information about solution and progress
	double LowerBound = -1.0;
	double UpperBound = 1e+20;
	double thresholdObjective = 1e+20;
	ArrayXXd BestSolution;
	bool UpdatedBound = false;

	double TimeLimit = 7200;
	double SolveTime = 0.0;
	int trustCount = -1;


	//Core point
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;	
	void UpdateCorePoint(const ArrayXXd& x_new);


	//Adding cuts
	inline void AddBendersCuts(const ArrayXXd& x_tilde, const ArrayXd& theta_tilde);
	inline void AddTrustCutsInner(IloModel& submodel, const ArrayXXd& x_tilde, double distance);
	inline void AddTrustCutsOuter(const ArrayXXd& x_tilde, double distance);
	inline void AddTabooCut(IloModel& submodel, const ArrayXXd& x_tilde);


	//Subproblems
	void IntegerCutCallback(ArrayXXd& x_tilde, ArrayXd& theta_tilde, double obj_tilde = -1, bool diversify = true, CallbackStatus status = CallbackStatus::FractionalRound);
	IloAlgorithm::Status SubproblemRestricted(const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj);
	IloAlgorithm::Status SubproblemDiversify(const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj);

	//Utility functions
	double CalculateObjective(const ArrayXXd& x_tilde);
	vector<VectorXd> CalculateItilde(const ArrayXXd& x_tilde);
	void UpdateSolution(const ArrayXXd& x_tilde, double obj);
	double GetUpperBound(IloCplex& cplex);
};

