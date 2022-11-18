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

#include "Model_Multicut.h"
#include "MulticutCallback.h"

#include "Greedy.h"


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

int main(int argc, char** argv) {

    bool verbose = true;
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
    std::string basefile = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC%d_Simple.json";

    json Stats;
    {
        std::ifstream f("C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\Overlap_Statistics.json", ifstream::in);
        f >> Stats;
        f.close();
    }
    map <string, map<string, vector<int>>> stats = Stats;




    //for (int test = 0; test < 10; test++) {
    //    time_t start;
    //    time(&start);
    //    Data data;

    //    string file = string_format(basefile, test);
    //    data.load(file, false);
    //    std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

    //    map <string, map<string, vector<int>>> stats;
    //    Model_Multicut mdl;
    //    mdl.SetData(data);

    //    mdl.threshold = 0.05;
    //    mdl.Solve(multicuts::Multi1SB1, useHeuristic::None, 0, true);
    //    cout << "M1SB1, threshold 0.05" << endl;
    //    cout << "Solution quality (data, optimal: 17989):" << data.SolutionQuality(mdl.Solution) << endl;
    //    cout << "Total time: " << mdl.TotalTime << endl;
    //    cout << "\n \n" << endl;


    //}






#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file;
    vector<std::string> datasets = {"Price"}; //Start = 16
    //vector<std::string> datasets = { "Distance" };
    //vector<std::string> datasets = {"HomeCharging"};

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

    for (auto dataset : datasets) {
        cout << "Dataset: " << dataset << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + dataset + "/";
        string prefix = "Final";
        std::string fp_times = resultspath + prefix + "_SolveTimes.json";
        std::string fp_gaps = resultspath + prefix + "_OptimalityGaps.json";
        std::string fp_objs = resultspath + prefix + "_ObjectiveValues.json";
        std::string fp_stats = resultspath + prefix + "_Statistics.json";

        

        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/MC%d.json";

        map<string, vector<double>> Results_times;
        {
            json Results;
            std::ifstream f(fp_times, ifstream::in);
            f >> Results;
            Results_times = Results;
            f.close();
        }

        map<string, vector<double>> Results_objs;
        {
            json Results;
            std::ifstream f(fp_objs, ifstream::in);
            f >> Results;
            Results_objs = Results;
            f.close();
        }

        map<string, vector<double>> Results_gaps;
        {
            json Results;
            std::ifstream f(fp_gaps, ifstream::in);
            f >> Results;
            f.close();
            Results_gaps = Results;
        }

        map<string, map<string, vector<int>>> Results_stats;     
        //{
        //    json Results;
        //    std::ifstream f(fp_stats, ifstream::in);
        //    f >> Results;
        //    f.close();
        //    Results_stats = Results;
        //}

        //vector<double> obj_a;
        //vector<double> time_a;
        //vector<double> gap_a;

        //vector<double> obj_b;
        //vector<double> time_b;
        //vector<double> gap_b;

        //vector<double> obj_c;
        //vector<double> time_c;
        //vector<double> gap_c;

        //vector<double> obj_d;
        //vector<double> time_d;
        //vector<double> gap_d;


        for (int test = 3; test < maxTest; test++) {
            cout << "Test: " << test << endl;
            time_t start;
            time(&start);
            Data data;
            file = string_format(basefile, test);
            data.load(file, verbose);
            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

            //{
            //    string label = "Multi1B1";
            //    Model_Multicut mdl;
            //    mdl.SetData(data);

            //    mdl.Solve(multicuts::Multi1B1, useHeuristic::None, 0);

            //    for (pair<string, int> res : mdl.stats) {
            //        string category = res.first;
            //        int value = res.second;
            //        vector<int>& ref = stats[label][category];
            //        ref.push_back(value);
            //    }

            //}

            {
                string label = "Multi1B1_WarmstartAndPostSimple";
                Model_Multicut mdl;
                mdl.SetData(data);

                mdl.Solve(multicuts::Multi1B1, useHeuristic::WarmstartAndPostSimple, 0);

                Results_objs[label].push_back(mdl.ObjectiveValue);
                Results_times[label].push_back(mdl.SolveTime);
                Results_gaps[label].push_back(mdl.OptimalityGap);
                for (pair<string, int> res : mdl.stats) {
                    string category = res.first;
                    int value = res.second;
                    Results_stats[label][category].push_back(value);
                }

            }


            {
                string label = "Multi1B1_WarmstartAndPostGreedy";
                Model_Multicut mdl;
                mdl.SetData(data);

                mdl.Solve(multicuts::Multi1B1, useHeuristic::WarmstartAndPostGreedy, 0);
                Results_objs[label].push_back(mdl.ObjectiveValue);
                Results_times[label].push_back(mdl.SolveTime);
                Results_gaps[label].push_back(mdl.OptimalityGap);
                for (pair<string, int> res : mdl.stats) {
                    string category = res.first;
                    int value = res.second;
                    Results_stats[label][category].push_back(value);
                }

            }

            {
                string label = "Multi3B2_WarmstartAndPostSimple";
                Model_Multicut mdl;
                mdl.SetData(data);

                mdl.Solve(multicuts::Multi3B2, useHeuristic::WarmstartAndPostSimple, 0);
                Results_objs[label].push_back(mdl.ObjectiveValue);
                Results_times[label].push_back(mdl.SolveTime);
                Results_gaps[label].push_back(mdl.OptimalityGap);
                for (pair<string, int> res : mdl.stats) {
                    string category = res.first;
                    int value = res.second;
                    Results_stats[label][category].push_back(value);
                }

            }

            {
                string label = "Multi3B2_WarmstartAndPostGreedy";
                Model_Multicut mdl;
                mdl.SetData(data);

                mdl.Solve(multicuts::Multi3B2, useHeuristic::WarmstartAndPostGreedy, 0);
                Results_objs[label].push_back(mdl.ObjectiveValue);
                Results_times[label].push_back(mdl.SolveTime);
                Results_gaps[label].push_back(mdl.OptimalityGap);
                for (pair<string, int> res : mdl.stats) {
                    string category = res.first;
                    int value = res.second;
                    Results_stats[label][category].push_back(value);
                }

            }

            //{
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    mdl.Solve(multicuts::Multi3B2, useHeuristic::WarmstartAndPostSimple, 0);
            //    obj_c.push_back(mdl.ObjectiveValue);
            //    time_c.push_back(mdl.SolveTime);
            //    gap_c.push_back(mdl.OptimalityGap);
            //    Results_times["Multi3B2_WarmstartAndPostSimple"] = time_c;
            //    Results_objs["Multi3B2_WarmstartAndPostSimple"] = obj_c;
            //    Results_gaps["Multi3B2_WarmstartAndPostSimple"] = gap_c;
            //}
            //{
            //    Model_Multicut mdl;
            //    mdl.SetData(data);
            //    mdl.Solve(multicuts::Multi3B2, useHeuristic::WarmstartAndPostGreedy, 0);
            //    obj_d.push_back(mdl.ObjectiveValue);
            //    time_d.push_back(mdl.SolveTime);
            //    gap_d.push_back(mdl.OptimalityGap);
            //    Results_times["Multi3B2_WarmstartAndPostGreedy"] = time_d;
            //    Results_objs["Multi3B2_WarmstartAndPostGreedy"] = obj_d;
            //    Results_gaps["Multi3B2_WarmstartAndPostGreedy"] = gap_d;
            //}



            cout << "\n" << endl;

            {
                std::ofstream f1(fp_times);
                f1 << std::setw(3) << (json) Results_times << std::endl;
                f1.close();
                std::ofstream f2(fp_objs);
                f2 << std::setw(3) << (json) Results_objs << std::endl;
                f2.close();
                std::ofstream f3(fp_gaps);
                f3 << std::setw(3) << (json) Results_gaps << std::endl;
                f3.close();
                std::ofstream f4(fp_stats);
                f4 << std::setw(3) << (json) Results_stats << std::endl;
                f4.close();
            }
        }
    }
#endif

	return 0;
}



