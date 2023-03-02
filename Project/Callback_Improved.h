#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data_Improved.h"


ILOSTLBEGIN
typedef IloArray<IloNumArray>    Num2D;
typedef IloArray<IloArray<IloNumArray>> Num3D;
typedef IloArray<IloNumVarArray>    NumVar2D;
typedef IloArray<NumVar2D>    NumVar3D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;
typedef IloArray<IloArray<IloArray<IloInt>>> Int3D;



class Callback_Improved : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	Callback_Improved();
	Callback_Improved(const Callback_Improved& tocopy);

	//Data members
	Data_Improved data;
	BoolVar2D x;
	IloNumVarArray theta;
	double incumbent = 0.0;
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;
	double obj_val = -1.0;
	short int duplicate_counter = 0;
	double dual_val = 1.7e307;

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	int M = data.params["M"];
	vector<int> Mj;


	//Functions for adding cuts
	void LazyCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);
	void UserCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);
	void AddCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde);



public:
	Callback_Improved(const Data_Improved& _data, const BoolVar2D& _x, const IloNumVarArray& _theta) :data(_data),  x(_x), theta(_theta) {
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
		for (int val : data.params["Mj"]) { Mj.push_back(val); }
	
	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = {};

	void UpdateCorePoint(ArrayXXd x_new);

	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;
};

