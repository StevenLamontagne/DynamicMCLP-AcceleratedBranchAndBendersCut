#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <chrono>

#include "Data.h"
#include "Utils.h"



/*
Class that handles the Benders cut generation as CPLEX callbacks. This method 
uses the optimality cuts as presented in Cordeau et al. (2019), without
any acceleration techniques. 
*/
class SingleCutBenders_Callback: public IloCplex::Callback::Function
{
private:
	SingleCutBenders_Callback();
	SingleCutBenders_Callback(const SingleCutBenders_Callback& tocopy);

	Data data;
	BoolVar2D x;
	IloNumVar theta;

	//Both of the following are used for looping detection
	int duplicate_counter = -1;
	double obj_val = 0.0;

	/*
	Function to run when the candidate solution is integer feasible.
	Calculates the I_tilde for this context.
	*/
	void LazyCutCallback(const IloCplex::Callback::Context& context);

	/*
	Function to run when the candidate solution is fractional.
	Calculates the I_tilde for this context.
	*/
	void UserCutCallback(const IloCplex::Callback::Context& context);

	/*
	Adds the Benders optimality cuts based on I_tilde. By default, this adds B1-type
	cuts, but this behaviour can be changed in the code. 
	*/
	void AddCuts(const IloCplex::Callback::Context& context, const vector<VectorXd>& I_tilde);


public:
	SingleCutBenders_Callback(const Data & _data, const BoolVar2D& _x, const IloNumVar & _theta):data(_data), x(_x), theta(_theta) {
	};

	/*
	Determines the context in which the callback has been activated (integer feasible or fractional solution), and invokes the 
	appropriate class functions. Also detects and fixes looping when it happens
	*/
	virtual void invoke(const IloCplex::Callback::Context& context);


	//Attributes for storing solve statistics. Since this class is initialised and 
	//deconstructed within the model class, they are not accessible generally. 
	map<string, int> stats = {};
	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;





};

