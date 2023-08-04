#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data.h"
#include "Utils.h"





/*
Class that handles the Benders cut generation as CPLEX callbacks. This method
uses Pareto-optimal cuts and multi-cut generation.
*/
class MultiCutBenders_Callback : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	MultiCutBenders_Callback();
	MultiCutBenders_Callback(const MultiCutBenders_Callback& tocopy);

	//Data members
	Data data;
	BoolVar2D x;
	IloNumVarArray theta;
	double incumbent = 0.0;

	//Used for looping detection
	double obj_val = -1.0;
	short int duplicate_counter = 0;

	//Core point information, for generating Pareto-optimal cuts
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;



	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;

	int M = data.params["M"];
	vector<int> Mj;


	/*
	Function to run when the candidate solution is integer feasible.
	Calculates the I_tilde for this context.
	*/
	void LazyCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);

	/*
	Function to run when the candidate solution is fractional.
	Calculates the I_tilde for this context.
	*/
	void UserCutCallback(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);

	/*
	Adds the Benders optimality cuts based on I_tilde. By default, this adds Pareto-optimal 
	B1-type cuts, but this behaviour can be changed in the code.
	*/
	void AddCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, const vector<VectorXd>& I_tilde);



public:
	MultiCutBenders_Callback(const Data& _data, const BoolVar2D& _x, const IloNumVarArray& _theta) :data(_data),  x(_x), theta(_theta) {
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
		for (int val : data.params["Mj"]) { Mj.push_back(val); }
	
	};

	/*
	Determines the context in which the callback has been activated (integer feasible or fractional solution), and invokes the
	appropriate class functions. Also detects and fixes looping when it happens
	*/
	virtual void invoke(const IloCplex::Callback::Context& context);

	/*
	Updates the core point and the coverage of the core point. This is only necessary if using Pareto-optimal cuts
	*/
	void UpdateCorePoint(ArrayXXd x_new);

	//Attributes for storing solve statistics. Since this class is initialised and 
	//deconstructed within the model class, they are not accessible generally. 
	map<string, int> stats = {};
	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;
};

