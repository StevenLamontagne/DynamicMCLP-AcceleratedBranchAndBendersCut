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

//#include "Data.h"
//#include "Model_Multicut.h"
//#include "MulticutCallback.h"
//#include "Greedy.h"
//
#include "Data_Improved.h"
#include "Greedy_Improved.h"
#include "Model_Improved.h"
//#include "Model_LocalBranching.h"

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
    final["Optimality gap"] = (double)temp.value("OptimalityGap (x100)", -1) / 100;
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
    std::string basefile = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\";

    //json Stats;
    //{
    //    std::ifstream f("C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\Overlap_Statistics.json", ifstream::in);
    //    f >> Stats;
    //    f.close();
    //}
    //map <string, map<string, vector<int>>> stats = Stats;

    json temp;
    vector<string> doubles = { "Solve time (sec)" , "Objective value" , "Optimality gap" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
    vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts" };
    vector<string> strings = { "Cplex status" };

    for (string key : doubles) { temp[key] = vector<double>(); }
    for (string key : ints) { temp[key] = vector<int>(); }
    for (string key : strings) { temp[key] = vector<string>(); }



    for (int test = 0; test < 1; test++) {
        time_t start;
        time(&start);

        Data_Improved data;
        data.load(basefile, string_format("MC%d_Simple.json", test), true);
        std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

        //Greedy_Improved G;
        //G.SetData(data);
        //G.Solve(false);

        //std::cout << "Greedy solving time: " << G.SolveTime << " seconds" << endl;
        //std::cout << "Greedy solution quality: " << G.SolutionQuality << endl;


        {
        string label = "Multi1B1_RootFractionalCallback";
        json params = {{ "verbose", true } };
        Model_Improved mdl;
        mdl.SetData(data);
        cout << "Method: " << label << endl;
        mdl.Solve(params);
        cout << "Solution quality (model): " << mdl.ObjectiveValue << endl;
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
        std::ofstream f1("Test.json");
        f1 << std::setw(3) << temp << std::endl;
        f1.close();


        //cout << "Solution quality (data): " << data.SolutionQuality(mdl.Solution) << endl;
        //std::ofstream f1("Solution.json");
        //json sol = { {"Exact", mdl.Solution}};
        //f1 << std::setw(3) << sol << std::endl;
        //f1.close();


        }





        //Data data;
        //data.load_compressed(basefile, string_format("MC%d_compressed.json", test), true);
        //std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

        //{
        //    string label = "Multi1B1_ImprovedTrust2";
        //    json params = { {"multicut", multicuts::Multi1B1}, {"heuristic", useHeuristic::Warmstart}, {"use_trust", true},  {"trust_threshold", 2}, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    cout << "Method: " << label << endl;
        //    mdl.Solve(params);

        //}

        //Model_BB mdl;
        //mdl.SetData(data);
        //mdl.Solve();
        //cout << "Objective value: " << mdl.ObjectiveValue << endl;
        //{
        //    json params = { {"multicut", multicuts::Multi1B1},{ "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    cout << "\n \n" << endl;
        //}

        //{
        //    json params = { {"multicut", multicuts::Multi3PO1}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    for (pair<string, int> res : mdl.stats) {
        //        string category = res.first;
        //        int value = res.second;
        //        cout << category + ": " << value << endl;
        //    }
        //    
        //}
        //{
        //    json params = { {"multicut", multicuts::Multi1B1}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    for (pair<string, int> res : mdl.stats) {
        //        string category = res.first;
        //        int value = res.second;
        //        cout << category + ": " << value << endl;
        //    }

        //}
        //{
        //    json params = { {"multicut", multicuts::SinglePO1}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    for (pair<string, int> res : mdl.stats) {
        //        string category = res.first;
        //        int value = res.second;
        //        cout << category + ": " << value << endl;
        //    }

        //}

        //{
        //    json params = { {"multicut", multicuts::Multi1PO1}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    for (pair<string, int> res : mdl.stats) {
        //        string category = res.first;
        //        int value = res.second;
        //        cout << category + ": " << value << endl;
        //    }

        //}
        //{
        //    json params = { {"multicut", multicuts::Multi1B1}, {"use_trust", true}, {"trust_threshold", 5.0 }, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    cout << "\n \n" << endl;
        //}

        //{
        //    json params = { {"multicut", multicuts::Multi1B1}, {"use_trust", true}, {"trust_threshold", 7.5 }, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    cout << "\n \n" << endl;
        //    cout << "Objective value: " << mdl.ObjectiveValue << endl;
        //    cout << "\n \n" << endl;
        //}

        //{
        //    json params = { {"multicut", multicuts::Multi1B1}, {"use_trust", true}, {"trust_threshold", 10.0 }, { "verbose", true } };
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve(params);
        //    cout << "\n \n" << endl;
        //}
        //{
        //    Model_Multicut mdl;
        //    mdl.SetData(data);
        //    mdl.Solve("", multicuts::Multi3B2, useHeuristic::Warmstart, 0, true);
        //    //cout << "Single B1" << endl;
        //    //cout << "Solution quality (data, optimal: 17989):" << data.SolutionQuality(mdl.Solution) << endl;
        //    //cout << "Total time: " << mdl.TotalTime << endl;
        //    cout << "\n \n" << endl;
        //}



    }






#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file;
    //vector<std::string> datasets = {"Price"}; //Start = 16
    //vector<std::string> datasets = { "Distance" };
    //vector<std::string> datasets = { "Simple" };
    //vector<std::string> datasets = {"LongSpan"};

    //for (auto dataset : datasets) {
    //    cout << "Dataset: " << dataset << endl;
    //    std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/MC%d.json";

    //    long int reduction_total = 0;
    //    long int max_total = 0;
    //    long int uncoverable_total = 0;
    //    long int precovered_total = 0;
    //    for (int test = 0; test < maxTest; test++) {
    //        long int reduction_test = 0;
    //        long int max_test = 0;
    //        long int uncoverable_test = 0;
    //        long int precovered_test = 0;
    //        file = string_format(basefile, test);
    //        data.load(file, false);
    //        for (int t = 0; t < data.T; t++) {
    //            for (int i = 0; i < data.N; i++) {
    //                for (int r1 = 0; r1 < data.R[i]; r1++) {
    //                    max_total += 1;
    //                    max_test += 1;
    //                    if (data.P[t][i][r1] == triplet::Uncoverable) {
    //                        uncoverable_test += 1;
    //                        uncoverable_total += 1;
    //                        break;
    //                    }
    //                    else if (data.P[t][i][r1] == triplet::Precovered) {
    //                        precovered_test += 1;
    //                        precovered_total += 1;
    //                        break;
    //                    }
    //                    else {
    //                        for (int r2 = data.R[i] - 1; r2 > r1; r2--) {
    //                            bool isMatch = true;
    //                            for (int j = 0; j < data.M; j++) {
    //                                for (int k = 0; k < data.Mj[j]; k++) {
    //                                    if (data.a[t][i][r1][j][k] != data.a[t][i][r2][j][k]) { isMatch = false; break; }
    //                                }
    //                                if (!isMatch) { break; }
    //                            }
    //                            if (isMatch) {
    //                                reduction_total += 1;
    //                                reduction_test += 1;
    //                                break;
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        cout << "Test " << test << ": " << endl;
    //        cout << "Clustering: " << 100 * (double)reduction_test / (double)max_test << " %" << endl;
    //        cout << "Uncoverable: " << 100 * (double)uncoverable_test / (double)max_test << " %" << endl;
    //        cout << "Precovered: " << 100 * (double)precovered_test / (double)max_test << " %" << endl;

    //    }
    //    cout << "Dataset average: " << endl;
    //    cout << "Clustering: " << 100 * (double)reduction_total / (double)max_total << " %" << endl;
    //    cout << "Uncoverable: " << 100 * (double)uncoverable_total / (double)max_total << " %" << endl;
    //    cout << "Precovered: " << 100 * (double)precovered_total / (double)max_total << " %" << endl;
    //    cout << "\n" << endl;
    //}
    string dataset = argv[1];
    //for (auto dataset : datasets) 
    {
        cout << "Dataset: " << dataset << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + dataset + "/";
        string prefix = "Eigen";
        std::string fp_times = resultspath + prefix + "_SolveTimes.json";
        std::string fp_gaps = resultspath + prefix + "_OptimalityGaps.json";
        std::string fp_objs = resultspath + prefix + "_ObjectiveValues.json";
        std::string fp_stats = resultspath + prefix + "_Statistics.json";

        
        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/";

        //map<string, vector<double>> Results_times;
        //{
        //    json Results;
        //    std::ifstream f(fp_times, ifstream::in);
        //    f >> Results;
        //    Results_times = Results;
        //    f.close();
        //}


        //map<string, vector<double>> Results_objs;
        //{
        //    json Results;
        //    std::ifstream f(fp_objs, ifstream::in);
        //    f >> Results;
        //    Results_objs = Results;
        //    f.close();
        //}


        //map<string, vector<double>> Results_gaps;
        //{
        //    json Results;
        //    std::ifstream f(fp_gaps, ifstream::in);
        //    f >> Results;
        //    f.close();
        //    Results_gaps = Results;
        //}


        //map<string, map<string, vector<int>>> Results_stats;     
        //{
        //    json Results;
        //    std::ifstream f(fp_stats, ifstream::in);
        //    f >> Results;
        //    f.close();
        //    Results_stats = Results;
        //}
        json Results_stats;
        //{
        //    json Results;
        //    std::ifstream f(fp_stats, ifstream::in);
        //    f >> Results;
        //    f.close();
        //    Results_stats = Results;
        //}

        vector<string> doubles = { "Solve time (sec)" , "Objective value" , "Optimality gap" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
        vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts" };
        vector<string> strings = { "Cplex status" };

        for (string key : doubles) { Results_stats[key] = vector<double>(); }
        for (string key : ints) { Results_stats[key] = vector<int>(); }
        for (string key : strings) { Results_stats[key] = vector<string>(); }



        for (int test = 0; test < maxTest; test++) {
        //for (int test = 0; test < 2; test++) {
            cout << "Test: " << test << endl;
            string folder = basefile + string_format("Test%d", test);

            time_t start;
            time(&start);
            Data_Improved data;
            data.load(basefile, string_format("MC%d_compressed.json", test), true);
            //Data data;
            //data.load_compressed(basefile, string_format("MC%d_compressed.json", test), true);

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
                string label = "Eigen_FractionalRoot";
                json params = { { "verbose", true } };
                Model_Improved mdl;
                mdl.SetData(data);
                cout << "Method: " << label << endl;
                mdl.Solve(params);
                json final = ConvertMap(mdl.stats);
                for (auto& el : final.items())
                {
                    Results_stats[el.key()].push_back(el.value());
                }
                cout << endl;

            }



            //{
            //    string label = "MultiB1";
            //    json params = { {"multicut", multicuts::Multi1B1}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "Multi1B1_ImprovedTrust2";
            //    json params = { {"multicut", multicuts::Multi1B1}, {"heuristic", useHeuristic::Warmstart}, {"use_trust", true},  {"trust_threshold", 2}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "SingleB1_Trust5";
            //    json params = { {"multicut", multicuts::SingleB1}, {"heuristic", useHeuristic::Warmstart}, {"use_trust", true},  {"trust_threshold", 5}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}


            //{
            //    string label = "SingleB1_Trust10";
            //    json params = { {"multicut", multicuts::SingleB1}, {"heuristic", useHeuristic::Warmstart}, {"use_trust", true},  {"trust_threshold", 10}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}


            //if (test > 12) {
            //    string label = "Multi1B1";
            //    json params = { {"multicut", multicuts::Multi1B1}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //if (test > 11) {
            //    string label = "Multi1B2";
            //    json params = { {"multicut", multicuts::Multi1B2}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "Multi2B1";
            //    json params = { {"multicut", multicuts::Multi2B1}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "Multi2B2";
            //    json params = { {"multicut", multicuts::Multi2B2}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "Multi3B1";
            //    json params = { {"multicut", multicuts::Multi3B1}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            //{
            //    string label = "Multi3B2";
            //    json params = { {"multicut", multicuts::Multi3B2}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}


            //{
            //    string label = "SingleB2_Warmstart";
            //    json params = { {"multicut", multicuts::SingleB2}, {"heuristic", useHeuristic::Warmstart}, { "verbose", true } };
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    cout << "Method: " << label << endl;
            //    mdl.Solve(params);
            //    Results_objs[label].push_back(mdl.ObjectiveValue);
            //    Results_times[label].push_back(mdl.SolveTime);
            //    Results_gaps[label].push_back(mdl.OptimalityGap);
            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        Results_stats[label][category].push_back(value);
            //    }

            //}

            cout << "\n" << endl;

            {
                //std::ofstream f1(fp_times);
                //f1 << std::setw(3) << (json) Results_times << std::endl;
                //f1.close();
                //std::ofstream f2(fp_objs);
                //f2 << std::setw(3) << (json) Results_objs << std::endl;
                //f2.close();
                //std::ofstream f3(fp_gaps);
                //f3 << std::setw(3) << (json) Results_gaps << std::endl;
                //f3.close();
                std::ofstream f4(fp_stats);
                f4 << std::setw(3) << (json) Results_stats << std::endl;
                f4.close();
            }
        }
    }
#endif

	return 0;
}



