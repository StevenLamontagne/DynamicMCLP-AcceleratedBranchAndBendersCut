#pragma once

#include <string>
#include <iostream>
#include <exception>
#include <memory>
#include <stdexcept>
#include <ctime>
#include <sys/stat.h>


#include "Librairies/nlohmann/json.hpp"
using json = nlohmann::json;


#include "Data.h"
#include "Greedy.h"
#include "BranchAndCut_Model.h"
#include "SingleCutBenders_Model.h"
#include "MultiCutBenders_Model.h"
#include "LocalBranching_Model.h"





//Simple map for converting the (integer-valued) results into more common forms
json ConvertMap(map<string, int> stats) {
    json temp = stats;
    json final;

    vector<string> StatusConversion = { "Unknown", "Feasible", "Optimal", "Infeasible", "Unbounded", "InfeasibleOrUnbounded", "Error", "Bounded"};
    final["Cplex status"] = StatusConversion[temp.value("CplexStatus", 0)];
    final["Solve time (sec)"] = (double)temp.value("SolveTime (x100)", -1) / 100;
    final["Objective value"] = (double)temp.value("ObjectiveValue (x100)", -1) / 100;
    final["Optimality gap (%)"] = (double)temp.value("OptimalityGap (x100)", -1) / 100;
    final["Number of nodes"] = temp.value("nNodes", -1);
    final["Number of lazy cuts"] = temp.value("nLazyCuts", -1);
    final["Average lazy cut time (sec)"] = (double)temp.value("LazyCutTime (x1000)", -1) / 1000;
    final["Number of user cuts"] = temp.value("nUserCuts", -1);
    final["Average user cut time (sec)"] = (double)temp.value("UserCutTime (x1000)", -1) / 1000;

    return final;
}

int main(int argc, char** argv) {

        std::string sharedFile = argv[1]; //Path to the file with dataset-level parameters (called 'Shared.json' in the repository)
        std::string instanceFile = argv[2]; //Path to the file with instance-level parameters (called 'MC_{instance}_compressed.json' in the repository)

        time_t start;
        time(&start);
        Data data;
        data.load(sharedFile, instanceFile, true);

        //If using the SepB method (which can only be used if there are two years), it suffices to comment out the following line
        //data.T = 2; 

        cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

        {
            string label = "Greedy";
            Greedy G;
            G.SetData(data);
            cout << "Method: " << label << endl;
            G.Solve(false, BUDGET_TYPE::Knapsack);

            std::cout << "Greedy solution quality: " << G.SolutionQuality << endl;
            std::cout << "Greedy solve time: " << G.SolveTime << endl << endl;
        }


        {
            string label = "B&C";
            json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            BranchAndCut_Model mdl;
            mdl.SetData(data);
            cout << "Method: " << label << endl;
            mdl.Solve(params);

            json final = ConvertMap(mdl.stats);
            for (auto& el : final.items())
            {
                std::cout << el.key() << ": " << el.value() << endl;
            }
            std::cout << endl;
        }


        {
            string label = "U-B&BC";
            json params = { { "verbose", true } };
            SingleCutBenders_Model mdl;
            mdl.SetData(data);
            cout << "Method: " << label << endl;
            mdl.Solve(params);

            json final = ConvertMap(mdl.stats);
            for (auto& el : final.items())
            {
                std::cout << el.key() << ": " << el.value() << endl;
            }
            std::cout << endl;
        }

        {
            string label = "A-B&BC";
            json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            MultiCutBenders_Model mdl;
            mdl.SetData(data);
            cout << "Method: " << label << endl;
            mdl.Solve(params);

            json final = ConvertMap(mdl.stats);
            for (auto& el : final.items())
            {
                std::cout << el.key() << ": " << el.value() << endl;
            }
            std::cout << endl;
        }

        {
            string label = "A-B&BC+LB-SubD-SepD";
            json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            LocalBranching_Model mdl;
            mdl.SetData(data);
            cout << "Method: " << label << endl;
            mdl.Solve(params);

            json final = ConvertMap(mdl.stats);
            for (auto& el : final.items())
            {
                std::cout << el.key() << ": " << el.value() << endl;
            }
            std::cout << endl;
        }



	return 0;
}



