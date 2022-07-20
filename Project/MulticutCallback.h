#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <chrono>

#include "Data.h"


ILOSTLBEGIN
typedef IloArray<IloNumArray>    Float2D;
typedef IloArray<IloNumVarArray>    NumVar2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;
typedef IloArray<IloArray<IloArray<IloNum>>> Num3D;
typedef IloArray<IloArray<IloArray<IloInt>>> Int3D;


enum class multicuts :short int {
	SingleB1,
	SingleB2,
	Multi1B1,
	Multi1B2,
	Multi2B1,
	Multi2B2,
	Multi3B1,
	Multi3B2
};

class MulticutCallback : public IloCplex::Callback::Function
{
private:
	MulticutCallback();
	MulticutCallback(const MulticutCallback& tocopy);

	Data data;
	BoolVar3D x;
	NumVar2D theta; // need array if multicut (need template?)
	multicuts cut_type = multicuts::SingleB2;


	//Define constants for easier reading and writing
	int& T = (data.T);
	int& M = (data.M);
	int& N = (data.N);
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;

	void LazyCutCallback(const IloCplex::Callback::Context& context);
	void UserCutCallback(const IloCplex::Callback::Context& context);
	void AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde);


public:
	MulticutCallback(const Data& _data, const BoolVar3D& _x, const NumVar2D& _theta, multicuts _cut_type) :data(_data), x(_x), theta(_theta), cut_type(_cut_type) {
	};

	virtual void invoke(const IloCplex::Callback::Context& context);







};


