#pragma once
#include <ilcplex/ilocplex.h>
#include <exception>
#include <ctime>

#include "Data.h"
#include "Utils.h"
#include "SingleCutBenders_Callback.h"


/*
Class for storing and solving the dynamic MCLP (EV charging station placement problem) using the 
U-B&BC method.
Solution value and objective value have dedicated attributes, whereas all other solve statistics (including
the solve time) are stored in the 'stats' attribute.
*/
class SingleCutBenders_Model
{
public:
	void SetData(const Data & newData) { data = newData; };
	void Solve(json params);

	//Options
	bool verbose = false;


	//Storing information about solution
	double ObjectiveValue = -1;
	ArrayXXd Solution;
	map<string, int> stats;

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

