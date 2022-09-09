#pragma once

#include <string>
#include <iostream>
#include <exception>
#include <memory>
#include <stdexcept>
#include <ctime>

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
	Data data;




#ifdef _WIN32
    std::string resultspath = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\Results\\";
    std::string file =  "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC0_Home.json";
    
    time_t start;
    time(&start);
    data.load(file, false);
    std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
    
    //long int reduction = 0;
    //long int total = 0;
    //for (int t = 0; t < data.T; t++) {
    //    for (int i = 0; i < data.N; i++) {
    //        for (int r1 = 0; r1 < data.R[i]; r1++) {
    //            total += 1;
    //            for (int r2 = data.R[i]-1; r2 > r1; r2--) {
    //                bool isMatch = true;
    //                for (int j = 0; j < data.M; j++) {
    //                    for (int k = 0; k < data.Mj[j]; k++) {
    //                        if (data.a[t][i][r1][j][k] != data.a[t][i][r2][j][k]) { isMatch = false; break; }
    //                    }
    //                    if (!isMatch) { break; }
    //                }
    //                if (isMatch) {
    //                    reduction += 1;
    //                    break;
    //                }
    //            }
    //        }
    //    }
    //}
    //cout << "Number of scenarios which can be reduced: " << reduction << endl;
    //cout << "This represets a reduction of " << (double)reduction / (double)total << " %" << endl;
    //cout << "\n" << endl;
    
    //vector<double> times_warmstart;
    ////vector<double> times_nowarmstart;
    Model_Multicut mdl;
    mdl.SetData(data);
    mdl.Solve(multicuts::Multi1B1, useHeuristic::None, 0);
    cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl.Solution) << endl;
    cout << "Total time: " << mdl.TotalTime << endl;
    cout << "Repair ratio: " << mdl.RepairRatio << endl;
    cout << "\n \n" << endl;

    //mdl.Solve(multicuts::Multi1B1, useHeuristic::WarmstartAndPostGreedy, 0);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Total time: " << mdl.TotalTime << endl;
    //cout << "Repair ratio: " << mdl.RepairRatio << endl;
    //cout << "\n \n" << endl;


    //mdl.Solve(multicuts::Multi1SB2, useHeuristic::None, 0);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Total time: " << mdl.TotalTime << endl;
    //cout << "\n \n" << endl;

    //mdl.Solve(multicuts::Multi1B1, useHeuristic::Warmstart, 50);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Total time: " << mdl.TotalTime << endl;
    //cout << "\n \n" << endl;

    //mdl.Solve(multicuts::Multi1B1, useHeuristic::None, 50);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Total time: " << mdl.TotalTime << endl;
    //cout << "\n \n" << endl;

    //mdl2.Solve(multicuts::SingleB2);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl2.Solution) << endl;
    //cout << "\n \n" << endl;

    //mdl2.Solve(multicuts::Multi1B1, useHeuristic::None);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl2.Solution) << endl;
    //cout << "\n \n" << endl;

    //mdl2.Solve(multicuts::Multi1B1, useHeuristic::Warmstart);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl2.Solution) << endl;
    //cout << "\n \n" << endl;

    //mdl2.Solve(multicuts::Multi1B1, useHeuristic::WarmstartAndPost);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl2.Solution) << endl;
    //cout << "\n \n" << endl;
    //mdl2.Solve(multicuts::Multi3B2);
    //cout << "Solution quality (data, optimal:18091.1):" << data.SolutionQuality(mdl2.Solution) << endl;
    //cout << "\n \n" << endl;
    //mdl2.Solve(multicuts::Multi1B1, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi1B1);
    //times_nowarmstart.push_back(mdl2.SolveTime);
    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;

    //mdl2.Solve(multicuts::Multi1B2, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi1B2);
    //times_nowarmstart.push_back(mdl2.SolveTime);

    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;

    //mdl2.Solve(multicuts::Multi2B1, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi2B1);
    //times_nowarmstart.push_back(mdl2.SolveTime);

    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;

    //mdl2.Solve(multicuts::Multi2B2, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi2B2);
    //times_nowarmstart.push_back(mdl2.SolveTime);

    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;

    //mdl2.Solve(multicuts::Multi3B1, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi3B1);
    //times_nowarmstart.push_back(mdl2.SolveTime);

    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;

    //mdl2.Solve(multicuts::Multi3B2, &mdl.Solution);
    ////std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    ////std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    ////std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //times_warmstart.push_back(mdl2.SolveTime);
    //mdl2.Solve(multicuts::Multi3B2);
    //times_nowarmstart.push_back(mdl2.SolveTime);


    //cout << "With warmstart: [" << endl;
    //for (auto i : times_warmstart) { cout << i<< ","; }
    //cout << "]" << endl;

    //cout << "Without warmstart: [" << endl;
    //for (auto i : times_nowarmstart) { cout << i << ","; }
    //cout << "]" << endl;
    //Results_times["CutTimes"] = *mdl2.CutTimes;
    //std::ofstream f2(resultspath + "CutTimes.json");
    //f2 << std::setw(3) << Results_times << std::endl;
    //f2.close();




#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file;
    //vector<std::string> datasets = {"HomeCharging"};
    vector<std::string> datasets = { "Price", "LongSpan" };
    //vector<std::string> datasets = {"Simple"};

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
        std::string fp_times = resultspath + "Final_SolveTimes.json";
        std::string fp_gaps = resultspath + "Final_OptimalityGaps.json";
        std::string fp_objs = resultspath + "Final_ObjectiveValues.json";

        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/MC%d.json";
        json Results_times;
        //{
        //    std::ifstream f(fp_times, ifstream::in);
        //    f >> Results_times;
        //    f.close();
        //}
        json Results_objs;
        //{
        //    std::ifstream f(fp_objs, ifstream::in);
        //    f >> Results_objs;
        //    f.close();
        //}
        json Results_gaps;
        //{
        //    std::ifstream f(fp_gaps, ifstream::in);
        //    f >> Results_gaps;
        //    f.close();
        //}
        vector<double> obj_a;
        vector<double> time_a;
        vector<double> gap_a;

        vector<double> obj_b;
        vector<double> time_b;
        vector<double> gap_b;

        vector<double> obj_c;
        vector<double> time_c;
        vector<double> gap_c;

        vector<double> obj_d;
        vector<double> time_d;
        vector<double> gap_d;



        //vector<double> obj_m3b2_10;
        //vector<double> time_m3b2_10;
        //vector<double> gap_m3b2_10;

        //vector<double> obj_m3b2_50;
        //vector<double> time_m3b2_50;
        //vector<double> gap_m3b2_50;

        //vector<double> obj_m3b2_ws10;
        //vector<double> time_m3b2_ws10;
        //vector<double> gap_m3b2_ws10;

        //vector<double> obj_m3b2_ws50;
        //vector<double> time_m3b2_ws50;
        //vector<double> gap_m3b2_ws50;


        for (int test = 0; test < maxTest; test++) {
            cout << "Test: " << test << endl;
            time_t start;
            time(&start);
            file = string_format(basefile, test);
            data.load(file, verbose);
            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
            {
                Model_Multicut mdl;
                mdl.SetData(data);
                mdl.Solve(multicuts::Multi1B1, useHeuristic::None, 0);
                obj_a.push_back(mdl.ObjectiveValue);
                time_a.push_back(mdl.SolveTime);
                gap_a.push_back(mdl.OptimalityGap);
                Results_times["Multi1B1"] = time_a;
                Results_objs["Multi1B1"] = obj_a;
                Results_gaps["Multi1B1"] = gap_a;
            }
            {
                Model_Multicut mdl;
                mdl.SetData(data);
                mdl.Solve(multicuts::Multi1B1, useHeuristic::Warmstart, 0);
                obj_b.push_back(mdl.ObjectiveValue);
                time_b.push_back(mdl.SolveTime);
                gap_b.push_back(mdl.OptimalityGap);
                Results_times["Multi1B1_Warmstart"] = time_b;
                Results_objs["Multi1B1_Warmstart"] = obj_b;
                Results_gaps["Multi1B1_Warmstart"] = gap_b;
            }
            {
                Model_Multicut mdl;
                mdl.SetData(data);
                mdl.Solve(multicuts::Multi3B2, useHeuristic::None, 0);
                obj_c.push_back(mdl.ObjectiveValue);
                time_c.push_back(mdl.SolveTime);
                gap_c.push_back(mdl.OptimalityGap);
                Results_times["Multi3B2"] = time_c;
                Results_objs["Multi3B2"] = obj_c;
                Results_gaps["Multi3B2"] = gap_c;
            }
            {
                Model_Multicut mdl;
                mdl.SetData(data);
                mdl.Solve(multicuts::Multi3B2, useHeuristic::Warmstart, 0);
                obj_d.push_back(mdl.ObjectiveValue);
                time_d.push_back(mdl.SolveTime);
                gap_d.push_back(mdl.OptimalityGap);
                Results_times["Multi3B2_Warmstart"] = time_d;
                Results_objs["Multi3B2_Warmstart"] = obj_d;
                Results_gaps["Multi3B2_Warmstart"] = gap_d;
            }

            //{
            //    mdl.Solve(multicuts::Multi3B2, useHeuristic::None, 10);
            //    obj_m3b2_10.push_back(mdl.ObjectiveValue);
            //    time_m3b2_10.push_back(mdl.SolveTime);
            //    gap_m3b2_10.push_back(mdl.OptimalityGap);
            //    Results_times["Multi3B2_GRASP10"] = time_m3b2_10;
            //    Results_objs["Multi3B2_GRASP10"] = obj_m3b2_10;
            //    Results_gaps["Multi3B2_GRASP10"] = gap_m3b2_10;
            //}
            //{
            //    mdl.Solve(multicuts::Multi3B2, useHeuristic::None, 50);
            //    obj_m3b2_50.push_back(mdl.ObjectiveValue);
            //    time_m3b2_50.push_back(mdl.SolveTime);
            //    gap_m3b2_50.push_back(mdl.OptimalityGap);
            //    Results_times["Multi3B2_GRASP50"] = time_m3b2_50;
            //    Results_objs["Multi3B2_GRASP50"] = obj_m3b2_50;
            //    Results_gaps["Multi3B2_GRASP50"] = gap_m3b2_50;
            //}
            //{
            //mdl.Solve(multicuts::Multi3B2, useHeuristic::Warmstart, 10);
            //obj_m3b2_ws10.push_back(mdl.ObjectiveValue);
            //time_m3b2_ws10.push_back(mdl.SolveTime);
            //gap_m3b2_ws10.push_back(mdl.OptimalityGap);
            //Results_times["Multi3B2_WS_GRASP10"] = time_m3b2_ws10;
            //Results_objs["Multi3B2_WS_GRASP10"] = obj_m3b2_ws10;
            //Results_gaps["Multi3B2_WS_GRASP10"] = gap_m3b2_ws10;
            //}
            //{
            //    mdl.Solve(multicuts::Multi3B2, useHeuristic::Warmstart, 50);
            //    obj_m3b2_ws50.push_back(mdl.ObjectiveValue);
            //    time_m3b2_ws50.push_back(mdl.SolveTime);
            //    gap_m3b2_ws50.push_back(mdl.OptimalityGap);
            //    Results_times["Multi3B2_WS_GRASP50"] = time_m3b2_ws50;
            //    Results_objs["Multi3B2_WS_GRASP50"] = obj_m3b2_ws50;
            //    Results_gaps["Multi3B2_WS_GRASP50"] = gap_m3b2_ws50;
            //}

            cout << "\n" << endl;

            {
                std::ofstream f1(fp_times);
                f1 << std::setw(3) << Results_times << std::endl;
                f1.close();
                std::ofstream f2(fp_objs);
                f2 << std::setw(3) << Results_objs << std::endl;
                f2.close();
                std::ofstream f3(fp_gaps);
                f3 << std::setw(3) << Results_gaps << std::endl;
                f3.close();
            }
        }
    }
#endif

	return 0;
}

