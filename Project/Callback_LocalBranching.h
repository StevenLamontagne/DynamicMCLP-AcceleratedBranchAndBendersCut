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
	double incumbent = 0.0;
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;
	double obj_val = -1.0;
	int obj_counter = 0;

	IloModel BaseModel;
	IloModel* MasterModel;

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	int M = data.params["M"];
	vector<int> Mj;


	//Functions for adding cuts
	void IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde, bool diversify = true);
	void FractionalCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde);

	void AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde);
	void AddTrustCuts(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, const double& minimum, const double& maximum);
	void AddTrustCuts(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, const double& maximum);

	double CalculateDistance(const ArrayXXd& x1, const ArrayXXd& x2);
	double CalculateDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x);
	double CalculateObjective(const ArrayXXd& x_tilde);
	vector<VectorXd> CalculateItilde(const ArrayXXd& x_tilde);



	IloAlgorithm::Status SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde,  ArrayXXd& sol, double& obj);
	IloAlgorithm::Status SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj);





public:
	Callback_LocalBranching(const Data_Improved& _data, IloModel* _model, const BoolVar2D& _x, const IloNumVarArray& _theta) :data(_data), MasterModel(_model), x(_x), theta(_theta) {
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
		for (int val : data.params["Mj"]) { Mj.push_back(val); }
		BaseModel = *_model;

	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = {};

	void UpdateCorePoint(ArrayXXd x_new);

	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;
};

