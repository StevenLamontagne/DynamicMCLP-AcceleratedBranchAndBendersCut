#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data_Improved.h"
#include "Callback_Improved.h"


ILOSTLBEGIN
typedef IloArray<IloNumArray>    Num2D;
typedef IloArray<IloArray<IloNumArray>> Num3D;
typedef IloArray<IloNumVarArray>    NumVar2D;
typedef IloArray<NumVar2D>    NumVar3D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;
typedef IloArray<IloArray<IloArray<IloInt>>> Int3D;


enum class CallbackStatus {
	Unknown,
	Optimal,
	Integer,
	FractionalFloor,
	FractionalRound
};

class Callback_LocalBranching : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	Callback_LocalBranching();
	Callback_LocalBranching(const Callback_LocalBranching& tocopy);

	//Data members
	Data_Improved data;
	BoolVar2D x;
	IloNumVarArray theta;
	IloCplex cplex;
	
	bool SkipFractional = false;

	
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


	//Core point calculation
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;

	
	

	

	vector<double> UserCutTimes;
};

