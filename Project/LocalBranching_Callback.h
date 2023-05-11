#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data.h"
#include "MultiCutBenders_Callback.h"




class LocalBranching_Callback : public IloCplex::Callback::Function
{
private:
	//Set to private to ban usage
	LocalBranching_Callback();
	LocalBranching_Callback(const LocalBranching_Callback& tocopy);

	//Data members
	Data data;
	IloModel MasterModel;
	BoolVar2D x;
	IloNumVarArray theta;
	IloNumVar theta_obj;
	BoolVar2D trust;


	int trustCount = -1;

	

	ArrayXXd core_point;
	vector<VectorXd> core_coverage;
	

	//Solution objects
	ArrayXXd BestSolution;
	double LowerBound = 0.0;
	bool UpdatedSolution = true;

	void UpdateSolution(const ArrayXXd& x_tilde, double obj)
	{
		BestSolution = x_tilde;
		LowerBound = obj;
		UpdatedSolution = true;
		cout << "Found solution of value: " << obj << endl;
	}

	
	

	//Define constants for easier reading and writing
	int& T = data.T;
	int& M_bar = data.M_bar;
	vector<int>& P = data.P;



	//Functions for adding cuts
	void IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde = -1, bool diversify = true);


	void AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);
	void AddTrustCutsInner(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, double distance);
	void AddTrustCutsMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, double distance);
	void AddTabooCut(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde);
	void AddTabooCutMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);


	bool WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2);
	bool WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x);

	double CalculateObjective(const ArrayXXd& x_tilde);
	vector<VectorXd> CalculateItilde(const ArrayXXd& x_tilde);

	void UpdateCorePoint(ArrayXXd x_new);

	IloAlgorithm::Status SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde,  ArrayXXd& sol, double& obj);
	IloAlgorithm::Status SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj);





public:
	LocalBranching_Callback(const Data& _data, const IloModel& _model, const BoolVar2D& _x, const IloNumVarArray& _theta, const IloNumVar& _theta_obj, const BoolVar2D& _trust) :data(_data), MasterModel(_model), x(_x), theta(_theta), theta_obj(_theta_obj), trust(_trust){
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
	};

	virtual void invoke(const IloCplex::Callback::Context& context);

	map<string, int> stats = { {"nRestricted",0}, {"nDiversified",0} };
	CALLBACK_STATUS status = CALLBACK_STATUS::Unknown;


	ArrayXXd GetSolution() { return BestSolution; }
	void SetSolution(ArrayXXd Solution, double Value) { BestSolution = Solution; LowerBound = Value; }
	double GetSolutionObjective() { return LowerBound; }

	
	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;
};

