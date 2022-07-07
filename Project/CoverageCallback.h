#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

#include "Data.h"
ILOSTLBEGIN
typedef IloArray<IloNumArray>    Float2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;
typedef IloArray<IloArray<IloArray<IloNum>>> Num3D;
typedef IloArray<IloArray<IloArray<IloInt>>> Int3D;


enum class cuts :char {
	SingleB0,
	SingleB1,
	SingleB2
};

class CoverageCallback: public IloCplex::Callback::Function
{
private:
	CoverageCallback();
	CoverageCallback(const CoverageCallback& tocopy);

	Data data;
	BoolVar3D x;
	IloNumVar theta;
	cuts cut_type = cuts::SingleB2;
	

	//Define constants for easier reading and writing
	int& T = (data.T);
	int& M = (data.M);
	int& N = (data.N);
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;


public:
	CoverageCallback(const Data & _data, const BoolVar3D& _x, const IloNumVar & _theta, cuts _cut_type):data(_data), x(_x), theta(_theta), cut_type(_cut_type) {
	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	void LazyCutCallback(const IloCplex::Callback::Context& context);
	void UserCutCallback(const IloCplex::Callback::Context& context);




};

