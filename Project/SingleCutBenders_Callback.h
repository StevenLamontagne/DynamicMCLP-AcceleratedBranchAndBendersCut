#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <chrono>

#include "Data.h"
#include "Utils.h"




class SingleCutBenders_Callback: public IloCplex::Callback::Function
{
private:
	SingleCutBenders_Callback();
	SingleCutBenders_Callback(const SingleCutBenders_Callback& tocopy);

	Data data;
	BoolVar2D x;
	IloNumVar theta;
	SINGLE_CUTS cut_type = SINGLE_CUTS::SingleB1;
	
	int duplicate_counter = -1;
	double obj_val = 0.0;

	void LazyCutCallback(const IloCplex::Callback::Context& context);
	void UserCutCallback(const IloCplex::Callback::Context& context);
	void AddCuts(const IloCplex::Callback::Context& context, const vector<VectorXd>& I_tilde);


public:
	SingleCutBenders_Callback(const Data & _data, const BoolVar2D& _x, const IloNumVar & _theta, SINGLE_CUTS _cut_type):data(_data), x(_x), theta(_theta), cut_type(_cut_type) {
	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = {};

	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;





};

