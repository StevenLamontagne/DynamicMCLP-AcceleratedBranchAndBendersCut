#pragma once

#include <string>
#include <iostream>
#include <exception>
#include <memory>
#include <stdexcept>
#include <ctime>
#include <sys/stat.h>


#ifdef _WIN32
#include <nlohmann/json.hpp>
#endif

#ifdef __linux__
#include "json.hpp"
//#include "jsonOld.hpp"
#endif


using json = nlohmann::json;


#include "Data.h"
#include "Greedy.h"
#include "BranchAndCut_Model.h"
#include "SingleCutBenders_Model.h"
#include "MultiCutBenders_Model.h"
#include "LocalBranching_Model.h"

template<typename ... Args>
std::string string_format(const std::string& format, Args ... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
    if (size <= 0) { throw std::runtime_error("Error during formatting."); }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

template <class _InIt1, class _InIt2>
double dot(_InIt1 const& vec1, _InIt2 const& vec2) { 
    if (vec1.size() != vec2.size()) { throw std::runtime_error("Supplied iterators are different sizes."); } 
    return std::inner_product(std::begin(vec1), std::end(vec1), std::begin(vec2), 0.0); 
}



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

    /*
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'v': verbose = true; break;
            default: std::cout << "Unrecognised argument \n";
            }
        }
        else
            file = argv[i];
    }
    */




#ifdef _WIN32
    std::string resultspath = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\Results\\";
    //std::string file = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC0_Home.json";

    //json Stats;
    //{
    //    std::ifstream f("C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\Overlap_Statistics.json", ifstream::in);
    //    f >> Stats;
    //    f.close();
    //}
    //map <string, map<string, vector<int>>> stats = Stats;

    json temp;
    vector<string> doubles = { "Solve time, LP (sec)" , "Solve time, MIP (sec)", "Objective value" , "Optimality gap (%)" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
    vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts", "Number of restricted subproblems", "Number of diversified subproblems" };
    vector<string> strings = { "Cplex status" };

    for (string key : doubles) { temp[key] = vector<double>(); }
    for (string key : ints) { temp[key] = vector<int>(); }
    for (string key : strings) { temp[key] = vector<string>(); }


    time_t start;
    time(&start);
    string path = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\";
    string shared = path + "Shared_Home.json";
    string instance = path + "MC0_Home_compressed.json";
    Data data;
    data.load(shared, instance, true);
    std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
    //std::cout << "Precovered: " << data.Precovered[0] + data.Precovered[1] + data.Precovered[2] + data.Precovered[3] << endl;

    {
        Greedy G;
        G.SetData(data);
        G.Solve(false, BUDGET_TYPE::Knapsack);

        std::cout << "Greedy solution quality: " << G.SolutionQuality << endl << endl;
    }


    {
        string label = "BranchAndCut";
        json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
        BranchAndCut_Model mdl;
        mdl.SetData(data);
        cout << "Method: " << label << endl;
        mdl.Solve(params);
        //for (pair<string, int> res : mdl.stats) {
        //    string category = res.first;
        //    int value = res.second;
        //    cout << category + ": " << value << endl;
        //}
        json final = ConvertMap(mdl.stats);
        for (auto& el : final.items())
        {
            temp[el.key()].push_back(el.value());
            std::cout << el.key() << ": " << el.value() << '\n';
        }
        std::cout << endl << endl;
    }

    {
        string label = "SingleCutBenders";
        json params = { { "verbose", true } };
        SingleCutBenders_Model mdl;
        mdl.SetData(data);
        cout << "Method: " << label << endl;
        mdl.Solve(params);
        //for (pair<string, int> res : mdl.stats) {
        //    string category = res.first;
        //    int value = res.second;
        //    cout << category + ": " << value << endl;
        //}
        json final = ConvertMap(mdl.stats);
        for (auto& el : final.items())
        {
            temp[el.key()].push_back(el.value());
            std::cout << el.key() << ": " << el.value() << '\n';
        }
        std::cout << endl << endl;
    }

    {
        string label = "MultiCutBenders";
        json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
        MultiCutBenders_Model mdl;
        mdl.SetData(data);
        cout << "Method: " << label << endl;
        mdl.Solve(params);
        //for (pair<string, int> res : mdl.stats) {
        //    string category = res.first;
        //    int value = res.second;
        //    cout << category + ": " << value << endl;
        //}
        json final = ConvertMap(mdl.stats);
        for (auto& el : final.items())
        {
            temp[el.key()].push_back(el.value());
            std::cout << el.key() << ": " << el.value() << '\n';
        }
        std::cout << endl << endl;
    }

    {
        string label = "LocalBranching";
        json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
        LocalBranching_Model mdl;
        mdl.SetData(data);
        cout << "Method: " << label << endl;
        mdl.Solve(params);
        //for (pair<string, int> res : mdl.stats) {
        //    string category = res.first;
        //    int value = res.second;
        //    cout << category + ": " << value << endl;
        //}
        json final = ConvertMap(mdl.stats);
        for (auto& el : final.items())
        {
            temp[el.key()].push_back(el.value());
            std::cout << el.key() << ": " << el.value() << '\n';
        }
        std::cout << endl << endl;
    }






#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file;
    string dataset = argv[1];
    {
        cout << "Dataset: " << dataset << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + dataset + "/";
        string prefix = "QuickTrust";
        std::string fp_times = resultspath + prefix + "_SolveTimes.json";
        std::string fp_gaps = resultspath + prefix + "_OptimalityGaps.json";
        std::string fp_objs = resultspath + prefix + "_ObjectiveValues.json";
        std::string fp_stats = resultspath + prefix + "_Statistics.json";

        
        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/";
        string sharedFile = basefile + "Shared.json";
        json Results_stats;
        {
            json Results;
            std::ifstream f(fp_stats, ifstream::in);
            f >> Results;
            f.close();
            Results_stats = Results;
        }

        vector<string> doubles = { "Solve time (sec)" , "Objective value" , "Optimality gap (%)" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
        vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts" };
        vector<string> strings = { "Cplex status" };
        vector<string> labels = { "Multi1PO1+Distance4" };


        for (string label : labels) {
            if (Results_stats.contains(label)){ continue; }
            Results_stats[label] = {};
            for (string key : doubles) { Results_stats[label][key] = vector<double>(); }
            for (string key : ints) { Results_stats[label][key] = vector<int>(); }
            for (string key : strings) { Results_stats[label][key] = vector<string>(); }
        }
        //Results_stats["Greedy"] = {};
        //Results_stats["Greedy"]["Solve time (sec)"] = vector<double>();
        //Results_stats["Greedy"]["Objective value"] = vector<double>();

        for (int test = 0; test < maxTest; test++) {
        //for (int test = 0; test < 2; test++) {
            cout << "Test: " << test << endl;
            string folder = basefile + string_format("Test%d", test);
            string instanceFile = basefile + string_format("MC%d_compressed.json", test);

            time_t start;
            time(&start);
            Data data;
            data.load(sharedFile, instanceFile, true);

            data.T = 2; //Remember to remove this later
            //data.load_fromUtilities(sharedFile, instanceFile, true, priceProfile::OnlyLevel3_Unperturbed);

            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

            //{
            //    string label = "B&C";
            //    Model_BB mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve();
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    Results_stats[label]["nNodes"].push_back(mdl.nNodes);
            //}
            {
                string label = "Multi1PO1+Distance4";
                if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
                if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
                json params = { { "verbose", true } };
                Model_Improved mdl;
                mdl.SetData(data);
                cout << "Method: " << label << endl;
                mdl.Solve(params);
                json final = ConvertMap(mdl.stats);
                for (auto& el : final.items())
                {
                    Results_stats[label][el.key()].push_back(el.value());
                }
                cout << endl;

            }



            //{
            //    string label = "Greedy";
            //    cout << "Method: " << label << endl;

            //    Greedy G;
            //    G.SetData(data);
            //    G.Solve(false, BUDGET_TYPE::Knapsack);

            //    Results_stats[label]["Solve time (sec)"].push_back(G.SolveTime);
            //    Results_stats[label]["Objective value"].push_back(G.SolutionQuality);
            //    cout << "Objective value: " << G.SolutionQuality << endl;
            //    cout << "Solve time: " << G.SolveTime << endl;
            //    cout << endl;

            //}

            //{
            //    string label = "BranchAndCut";
            //    if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
            //    if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
            //    json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            //    BranchAndCut_Model mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    json final = ConvertMap(mdl.stats);
            //    //for (auto& el : Results_stats[label].items()) {
            //    //    if (!final.contains(el.key())) { Results_stats[label].erase(el.key()); }
            //    //}
            //    for (auto& el : final.items())
            //    {
            //        Results_stats[label][el.key()].push_back(el.value());
            //    }
            //    cout << endl;

            //}

            //{
            //    string label = "SingleCutBenders";
            //    if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
            //    if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
            //    json params = { { "verbose", true } };
            //    SingleCutBenders_Model mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    json final = ConvertMap(mdl.stats);
            //    for (auto& el : Results_stats[label].items()) {
            //        if (!final.contains(el.key())) { Results_stats[label].erase(el.key()); }
            //    }
            //    for (auto& el : final.items())
            //    {
            //        Results_stats[label][el.key()].push_back(el.value());
            //    }
            //    cout << endl;

            //}


            //{
            //    string label = "MultiCutBenders";
            //    if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
            //    if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
            //    json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            //    MultiCutBenders_Model mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    json final = ConvertMap(mdl.stats);
            //    for (auto& el : Results_stats[label].items()) {
            //        if (!final.contains(el.key())) { Results_stats[label].erase(el.key()); }
            //    }
            //    for (auto& el : final.items())
            //    {
            //        Results_stats[label][el.key()].push_back(el.value());
            //    }
            //    cout << endl;

            //}

            {
                string label = "LocalBranching_Cutting";
                if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
                if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
                json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
                LocalBranching_Model mdl;
                mdl.SetData(data);
                cout << "Method: " << label << endl;
                mdl.Solve(params);
                json final = ConvertMap(mdl.stats);
                for (auto& el : Results_stats[label].items()) {
                    if (!final.contains(el.key())) { Results_stats[label].erase(el.key()); }
                }
                for (auto& el : final.items())
                {
                    Results_stats[label][el.key()].push_back(el.value());
                }
                cout << endl;
            }


            //{
            //    Data data;
            //    //data.load(sharedFile, instanceFile, true);
            //    data.load_withLevel2(sharedFile, instanceFile, true, priceProfile::MixedLevel2andLevel3_Perturbed);

            //    cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

            //    {
            //        string label = "Perturbed_Greedy";
            //        Greedy G;
            //        G.SetData(data);
            //        G.Solve(false, BUDGET_TYPE::Knapsack);

            //        Results_stats[label]["Solve time (sec)"].push_back(G.SolveTime);
            //        Results_stats[label]["Objective value"].push_back(G.SolutionQuality);
            //    }
            //    {
            //        string label = "Perturbed_LocalBranching";
            //        if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
            //        if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
            //        json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack}, {"nGRASP", 0} };
            //        LocalBranching_Model mdl;
            //        mdl.SetData(data);
            //        cout << "Method: " << label << endl;
            //        mdl.Solve(params);
            //        json final = ConvertMap(mdl.stats);
            //        for (auto& el : final.items())
            //        {
            //            Results_stats[label][el.key()].push_back(el.value());
            //        }
            //        cout << endl;

            //    }
            //}

            cout << "\n" << endl;

            {
                std::ofstream f4(fp_stats);
                f4 << std::setw(3) << (json) Results_stats << std::endl;
                f4.close();
            }
        }
    }
#endif

	return 0;
}



