#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <map>
#include <chrono>


#include "Data.h"
#include "MultiCutBenders_Callback.h"



/*
Class that handles the Benders cut generation as CPLEX callbacks. This method
uses Pareto-optimal cuts and multi-cut generation, as well as implementing 
local branching. Uses SepD separation method
*/
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

	//Tracks what index to use for the binary auxiliary variables
	//in the separation
	int trustCount = -1;

	
	//Core point information, for generating Pareto-optimal cuts
	ArrayXXd core_point;
	vector<VectorXd> core_coverage;
	

	//Solution objects
	ArrayXXd BestSolution; //Due to the trust cuts, CPLEX is likely to not have the optimal solution array. It is important to keep it stored here.
	double LowerBound = 0.0;
	bool UpdatedSolution = true;

	/*
	Updates the incumbent and sets a flag to send a post-heuristic to CPLEX
	*/
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



	/*
	Function to run when the candidate solution is integer feasible.
	If the true objective value is known, it is passed as 'obj_tilde'. Otherwise it is calculated at the start of 
	the function. 
	The 'diversify' parameter controls whether or not to attempt a diversified subproblem around 'x_tilde' if 
	the restricted subproblem fails. 
	*/
	void IntegerCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde, double obj_tilde = -1, bool diversify = true);


	/*
	Add optimality cuts to the main problem
	*/
	void AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);

	/*
	Add cuts to model 'model' imposing that the distance from 'x_tilde' must be <= 'distance'.
	Used in the restricted and diversified subproblems.
	*/
	void AddTrustCutsInner(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde, double distance);

	/*
	Separates the subdomain of distance 'distance' around 'x_tilde' from the main problem. This is done via
	big-M constraints and binary variables. If the number of constraints has exceeded the maximum, these cuts are
	automatically replaced by taboo/no-good cuts, and a message is displayed. 
	*/
	void AddTrustCutsMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde, double distance);

	/*
	Separates the subdomain of distance 1 around 'x_tilde' from model 'model'. 
	Used in the diversified subproblems.
	*/
	void AddTabooCut(IloEnv& env, IloModel& model, const ArrayXXd& x_tilde);

	/*
	Separates the subdomain of distance 1 around 'x_tilde' from the main problem.
	*/
	void AddTabooCutMaster(IloEnv& env, const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde);

	/*
	Checks if the the distance between 'x1' and 'x2' (using our metric) is <= 2.
	This function is used when there is no guarantee that x1 >= x2 pointwise.
	*/
	bool WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2);

	/*
	Checks if the the distance between 'x_upper' and 'x' (using our metric) is <= 2.
	This function is used when x_upper >= x pointwise.
	*/
	bool WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x);

	/*
	Calculates the true objective value for 'x_tilde', since CPLEX may
	(unknowingly) be lying to us.
	*/
	double CalculateObjective(const ArrayXXd& x_tilde);

	/*
	Calculate the associated value of I_tilde given 'x_tilde'.
	*/
	vector<VectorXd> CalculateItilde(const ArrayXXd& x_tilde);

	/*
	Updates the core point and the coverage of the core point. This is only necessary if using Pareto-optimal cuts
	*/
	void UpdateCorePoint(ArrayXXd x_new);


	/*
	Solve the subproblem of distance 2 around 'x_tilde'. Returns the true objective value in 'obj', which is 
	used for the integer callback.
	*/
	IloAlgorithm::Status SubproblemRestricted(IloEnv& env, const ArrayXXd& x_tilde,  ArrayXXd& sol, double& obj);

	/*
	Solve the subproblem of distance 1 to 4 around 'x_tilde'. Returns the true objective value in 'obj', which is
	used for the integer callback.
	*/
	IloAlgorithm::Status SubproblemDiversify(IloEnv& env, const ArrayXXd& x_tilde, ArrayXXd& sol, double& obj);





public:
	LocalBranching_Callback(const Data& _data, const IloModel& _model, const BoolVar2D& _x, const IloNumVarArray& _theta, const IloNumVar& _theta_obj, const BoolVar2D& _trust) :data(_data), MasterModel(_model), x(_x), theta(_theta), theta_obj(_theta_obj), trust(_trust){
		core_point = ArrayXXd::Constant(T, M_bar, 0.0);
		for (int t = 0; t < T; t++) {
			VectorXd cover = VectorXd::Constant(P[t], 0.0);
			core_coverage.push_back(cover);
		}
	};

	/*
	Determines the context in which the callback has been activated (integer feasible or fractional solution), and invokes the
	appropriate class functions. Also updates the lower bound if a better incumbent was found.
	*/
	virtual void invoke(const IloCplex::Callback::Context& context);




	ArrayXXd GetSolution() { return BestSolution; }
	void SetSolution(ArrayXXd Solution, double Value) { BestSolution = Solution; LowerBound = Value; }
	double GetSolutionObjective() { return LowerBound; }

	//Attributes for storing solve statistics. Since this class is initialised and 
	//deconstructed within the model class, they are not accessible generally. 
	map<string, int> stats = { {"nRestricted",0}, {"nDiversified",0} };
	vector<double> LazyCutTimes;
	vector<double> UserCutTimes;
};

