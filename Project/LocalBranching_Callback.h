#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data.h"
#include "Utils.h"
#include "MultiCutBenders_Callback.h"



<<<<<<< Updated upstream:Project/Callback_LocalBranching.h

enum class CallbackStatus {
	Unknown,
	Optimal,
	Integer,
	FractionalFloor,
	FractionalRound
};
=======
>>>>>>> Stashed changes:Project/LocalBranching_Callback.h

class LocalBranching_Callback : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	LocalBranching_Callback();
	LocalBranching_Callback(const LocalBranching_Callback& tocopy);

	//Data members
<<<<<<< Updated upstream:Project/Callback_LocalBranching.h
	Data_Improved data;
	BoolVar2D x;
	IloNumVarArray theta;
	IloCplex cplex;
	
	bool SkipFractional = false;
=======
	Data data;
	IloModel MasterModel;
	BoolVar2D x;
	IloNumVarArray theta;
	IloNumVar theta_obj;
	BoolVar2D trust;


	int trustCount = -1;
>>>>>>> Stashed changes:Project/LocalBranching_Callback.h

	
	//Define constants for easier reading and writing
	int T = data.T;
	int M_bar = data.M_bar;
	vector<int>& P = data.P;



	//Functions for adding cuts
	void FractionalCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde);
	void AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);

	bool WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2);
	bool WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x);

	double CalculateObjective(const ArrayXXd& x_tilde);
	vector<VectorXd> CalculateItilde(const ArrayXXd& x_tilde);


public:
<<<<<<< Updated upstream:Project/Callback_LocalBranching.h
	Callback_LocalBranching(const Data_Improved& _data, const IloCplex& _cplex, const BoolVar2D& _x, const IloNumVarArray& _theta) :data(_data),  cplex(_cplex), x(_x), theta(_theta) {};

	virtual void invoke(const IloCplex::Callback::Context& context);

	//Solution objects
	map<string, int> stats = {};
	CallbackStatus status = CallbackStatus::Unknown;
	ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
	ArrayXd theta_tilde = ArrayXd::Zero(T);
	double thresholdObjective = 0.0;
	double LowerBound = 0.0;
	double UpperBound = 1e+20;
=======
	LocalBranching_Callback(const Data& _data, const IloModel& _model, const BoolVar2D& _x, const IloNumVarArray& _theta, const IloNumVar& _theta_obj, const BoolVar2D& _trust) :data(_data), MasterModel(_model), x(_x), theta(_theta), theta_obj(_theta_obj), trust(_trust){
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = { {"nRestricted",0}, {"nDiversified",0} };
	CALLBACK_STATUS status = CALLBACK_STATUS::Unknown;
>>>>>>> Stashed changes:Project/LocalBranching_Callback.h


	//Core point calculation
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;

	
	

	

	vector<double> UserCutTimes;
};

