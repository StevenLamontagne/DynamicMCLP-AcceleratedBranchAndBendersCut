#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>
#include <chrono>
#include <map>

#include "Data_Strengthened.h"


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
	Multi3B2,
	Multi1SB1,
	Multi1SB2,
	Multi3SB1,
	Multi3SB2

};

enum class useHeuristic :char {
	None,
	Warmstart,
	WarmstartAndPostGreedy,
	WarmstartAndPostSimple
};

class MulticutCallback_Strengthened : public IloCplex::Callback::Function
{
private:
	MulticutCallback_Strengthened();
	MulticutCallback_Strengthened(const MulticutCallback_Strengthened& tocopy);

	Data_Strengthened data;
	BoolVar3D z;
	BoolVar2D y;
	NumVar2D theta; // need array if multicut (need template?)
	multicuts cut_type = multicuts::SingleB2;
	useHeuristic heuristic = useHeuristic::Warmstart;
	double incumbent = 0.0;
	IloArray<IloExpr> BudgetConstraints;


	//Define constants for easier reading and writing
	int& T = (data.T);
	int& M = (data.M);
	int& N = (data.N);
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;

	void LazyCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde);
	void UserCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde);
	void AddCuts(const IloCplex::Callback::Context& context, const Num3D& I_tilde);


public:
	MulticutCallback_Strengthened(const Data_Strengthened& _data, const BoolVar3D& _z, const BoolVar2D& _y, const NumVar2D& _theta, multicuts _cut_type, useHeuristic _heuristic) :data(_data), z(_z), y(_y), theta(_theta), cut_type(_cut_type), heuristic(_heuristic) {

	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	int repairCounter = 0;
	int totalCounter = 0;



};


