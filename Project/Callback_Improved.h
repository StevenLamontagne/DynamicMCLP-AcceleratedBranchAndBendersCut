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
	vector<vector<double>> core_point;
	double obj_val = -1.0;
	int obj_counter = 0;

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	int M = data.params["M"];
	vector<int> Mj;


	//Functions for adding cuts
	void LazyCutCallback(const IloCplex::Callback::Context& context, const Num2D& x_tilde);
	void UserCutCallback(const IloCplex::Callback::Context& context, const Num2D& x_tilde);
	void AddCuts(const IloCplex::Callback::Context& context, const Num2D& x_tilde, const vector<VectorXd>& I_tilde);



public:
	Callback_Improved(const Data_Improved& _data, const BoolVar2D& _x, const IloNumVarArray& _theta) :data(_data), x(_x), theta(_theta) {
		for (int t = 0; t < T; t++) {
			vector<double> cp1;
			for (int j_bar = 0; j_bar < M_bar; j_bar++) {
				cp1.push_back(0.0);
			}
			core_point.push_back(cp1);
		}
		for (int val : data.params["Mj"]) { Mj.push_back(val); }

	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = {
		{"nOverlapCuts", 0},
		{"nPostRepairs", 0},
		{"nPostFills", 0},
		{"nPostInfeasible", 0},
		{"nLazyCuts", 0},
		{"nUserCuts", 0}
	};
};

