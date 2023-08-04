#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>


#include "Data.h"
#include "Utils.h"
#include "MultiCutBenders_Callback.h"
#include "LocalBranching_Callback.h"
#include "LocalBranching_TwoYearCallback.h"
#include "Greedy.h"



/*
Class for storing and solving the dynamic MCLP (EV charging station placement problem) using the
A-B&BC+LB methods.
Solution value and objective value have dedicated attributes, whereas all other solve statistics (including
the solve time) are stored in the 'stats' attribute.

'budgetType' is an option for using either knapsack-style or cardinality-style budget constraints. This was not
in the final paper, but was left in case it could be useful reference.
*/
class LocalBranching_Model
{
public:
	void SetData(const Data& newData);
	void Solve(json params);



	//Storing information about solution
	double ObjectiveValue = -1;
	ArrayXXd Solution;
	ArrayXd ValueByYear;
	map<string, int> stats;

	//Options
	bool verbose = false;
	//The array of binary variables for the disjuctive separation method must be preallocated
	//This parameter determines the maximum number of separations that can be performed. 
	//Any in excess of this value will use no-good cuts instead
	int nTrust = 1024; 
	BUDGET_TYPE budgetType = BUDGET_TYPE::Knapsack;

	

	BoolVar2D x;
	IloNumVarArray theta;


private:
	Data data;


	void GetSolution(IloCplex& cplex, BoolVar2D& x) {
		for (int t = 0; t < data.T; t++) {
			for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
				Solution(t, j_bar) = cplex.getValue(x[t][j_bar]);
			}
		}
	};

	void ConvertSolution() {
		ArrayXXd temp = ArrayXXd::Constant(data.T, data.params["M"], 0.0);
		for (int t = 0; t < data.T; t++) {
			for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
				int j = data.params["station_coord"][j_bar][0];
				int k = data.params["station_coord"][j_bar][1];
				if ((Solution(t, j_bar) > 1 - EPS) && (temp(t, j) < k)) { temp(t, j) = k; }
			}
		}

		Solution = temp;
	};

};

