#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
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
	Multi3B2,
	Multi1SB1,
	Multi3SB2

};

enum class useHeuristic :char {
	None,
	Warmstart,
	WarmstartAndPostGreedy,
	WarmstartAndPostSimple
};

class MulticutCallback : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	MulticutCallback();
	MulticutCallback(const MulticutCallback& tocopy);

	//Data members
	Data data;
	BoolVar3D x;
	NumVar2D theta; 
	multicuts cut_type = multicuts::Multi1B1;
	useHeuristic heuristic = useHeuristic::Warmstart;
	double incumbent = 0.0;

	//Integer overflow is intended here
	unsigned char PHcounter = 1;


	//Define constants for easier reading and writing
	int& T = (data.T);
	int& M = (data.M);
	int& N = (data.N);
	vector<int>& Mj = data.Mj;
	vector<int>& R = data.R;

	//Functions for adding cuts
	void LazyCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde);
	void UserCutCallback(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde);
	void AddCuts(const IloCplex::Callback::Context& context, const IloArray<IloArray<IloNumArray>>& x_tilde, const Num3D& I_tilde);

	//Functions for heuristics
	pair<int, int> argmax(map<pair<int, int>, double> vec);
	vector<map<int, vector<int>>> GetFractional(const IloArray<IloArray<IloNumArray>>& x_tilde);
	bool SimpleRepair(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage);
	bool GreedyRepair(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage);
	bool GreedyFill(IloEnv& env, IloArray<IloArray<IloNumArray>>& x_tilde, vector<double>& Budget, IloArray<IloArray<IloBoolArray>>& coverage);
	void AddPostHeuristic(const IloCplex::Callback::Context& context, IloArray<IloArray<IloNumArray>>& x_tilde);


public:
	MulticutCallback(const Data& _data, const BoolVar3D& _x, const NumVar2D& _theta,multicuts _cut_type, useHeuristic _heuristic) :data(_data), x(_x), theta(_theta), cut_type(_cut_type), heuristic(_heuristic) {

	};
	virtual void invoke(const IloCplex::Callback::Context& context);
	double threshold = 0.1;
	map<string, int> stats = {
		{"nOverlapCuts", 0},
		{"nPostRepairs", 0},
		{"nPostFills", 0},
		{"nPostInfeasible", 0},
		{"nLazyCuts", 0},
		{"nUserCuts", 0}
	};





};


