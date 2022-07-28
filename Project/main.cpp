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
#endif


using json = nlohmann::json;

#include "Data.h"

#include "Model_BB.h"

//#include "Model_BendersCordeau.h"
//#include "CoverageCallback.h"

#include "Model_Multicut.h"
#include "MulticutCallback.h"

//#include "Greedy.h"


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
    std::string file =  "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC1.json";
    
    time_t start;
    time(&start);
    data.load(file, verbose);
    std::cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
    vector<double> times_warmstart;
    vector<double> times_nowarmstart;


    //Greedy mdl;
    //mdl.SetData(data);
    //mdl.Solve();
    //cout << "Solution value (Greedy): " << mdl.SolutionQuality << endl;
    //cout << "Solution value (Data): " << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Solve time (seconds): " << mdl.SolveTime << endl;

    //std::cout << "\n" << endl;
    //std::cout << "\n" << endl;


    //Model_Multicut mdl2;
    //mdl2.SetData(data);

    /*
    vector<vector<double>> testSol = { {0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                                {0.2, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                                {0.6, 0.5, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
    */

    vector<vector<double>> testSol = { {0.1, 0.1, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                        {0.1, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                        {0.2, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
                        {0.6, 0.5, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
    vector<double> Budget = data.RemainingBudget(testSol);
  
    vector<vector<vector<bool>>> coverage;
    //Initialise coverage of triplets
    for (int t = 0; t < data.T; t++) {
        vector<vector<bool>> cover1;
        for (int i = 0; i < data.N; i++) {
            vector<bool> cover2;
            double weight = (double)data.params["Ni"][t][i] / (double)data.R[i];
            for (int r = 0; r < data.R[i]; r++) {
                bool val = 0;
                switch (data.P[t][i][r])
                {
                case triplet::Uncoverable:
                    val = 1; //Set to 1 to skip trying to cover later
                    break;
                case triplet::Precovered:
                    val = 1;
                    break;
                default:
                    vector<pair<int, int>> cover = data.cover[t][i][r];
                    for (int jBar = 0; jBar < cover.size(); jBar++) {
                        int j = cover[jBar].first;
                        int k0 = cover[jBar].second;
                        //NOTE: This uses the value of the original solution, but rounded down. This is deliberate, since we know that
                        //the new solution will be guaranteed to have at least these values
                        if (testSol[t][j] >= k0) {
                            val = 1;
                            break;
                        }
                    }
                    break;
                }
                cover2.push_back(val);
            }
            cover1.push_back(cover2);
        }
        coverage.push_back(cover1);
    }

    mdl2.GreedyRepair(testSol, Budget, coverage);
    cout << "Repaired solution" << endl;
    for (int t = 0; t < testSol.size(); t++) {
        cout << "Year " << t << endl;
        for (int j = 0; j < testSol[t].size(); j++) {
            if (testSol[t][j] > 0) {
                cout << "Station " << j << " :" << testSol[t][j] << endl;
            }
        }
    }

    cout << "\n" << endl;

    mdl2.GreedyFill(testSol, Budget, coverage);
    cout << "Filled solution" << endl;
    for (int t = 0; t < testSol.size(); t++) {
        cout << "Year " << t << endl;
        for (int j = 0; j < testSol[t].size(); j++) {
            if (testSol[t][j] > 0) {
                cout << "Station " << j << " :" << testSol[t][j] << endl;
            }
        }
    }

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
    vector<std::string> datasets = {"HomeCharging", "Distance" };
    //vector<std::string> datasets = { "Distance" };


    for (long unsigned int dBar = 0; dBar < datasets.size(); dBar++) {
        cout << "Dataset: " << datasets[dBar] << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + datasets[dBar] + "/";
        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + datasets[dBar] + "/MaximumCover/MC%d.json";
        json Results_times;
        {
            std::ifstream f(resultspath + "SolveTimes.json", ifstream::in);
            f >> Results_times;
            f.close();
        }
        json Results_objs;
        {
            std::ifstream f(resultspath + "ObjectiveValues.json", ifstream::in);
            f >> Results_objs;
            f.close();
        }
        json Results_gaps;
        {
            std::ifstream f(resultspath + "OptimalityGaps.json", ifstream::in);
            f >> Results_gaps;
            f.close();
        }


        vector<double> obj_m1b1;
        vector<double> time_m1b1;
        vector<double> gap_m1b1;

        vector<double> obj_m1b2;
        vector<double> time_m1b2;
        vector<double> gap_m1b2;

        vector<double> obj_m2b1;
        vector<double> time_m2b1;
        vector<double> gap_m2b1;

        vector<double> obj_m2b2;
        vector<double> time_m2b2;
        vector<double> gap_m2b2;

        vector<double> obj_m3b1;
        vector<double> time_m3b1;
        vector<double> gap_m3b1;

        vector<double> obj_m3b2;
        vector<double> time_m3b2;
        vector<double> gap_m3b2;

        for (int test = 0; test < maxTest; test++) {
            cout << "Test: " << test << endl;
            time_t start;
            time(&start);
            file = string_format(basefile, test);
            data.load(file, verbose);
            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
            Model_Multicut mdl;
            mdl.SetData(data);
            {
                mdl.Solve(multicuts::Multi1B1);
                obj_m1b1.push_back(mdl.ObjectiveValue);
                time_m1b1.push_back(mdl.SolveTime);
                gap_m1b1.push_back(mdl.OptimalityGap);
                Results_times["Multi1B1"] = time_m1b1;
                Results_objs["Multi1B1"] = obj_m1b1;
                Results_gaps["Multi1B1"] = gap_m1b1;
            }
            {
                mdl.Solve(multicuts::Multi1B2);
                obj_m1b2.push_back(mdl.ObjectiveValue);
                time_m1b2.push_back(mdl.SolveTime);
                gap_m1b2.push_back(mdl.OptimalityGap);
                Results_times["Multi1B2"] = time_m1b2;
                Results_objs["Multi1B2"] = obj_m1b2;
                Results_gaps["Multi1B2"] = gap_m1b2;
            }
            {
                mdl.Solve(multicuts::Multi2B1);
                obj_m2b1.push_back(mdl.ObjectiveValue);
                time_m2b1.push_back(mdl.SolveTime);
                gap_m2b1.push_back(mdl.OptimalityGap);
                Results_times["Multi2B1"] = time_m2b1;
                Results_objs["Multi2B1"] = obj_m2b1;
                Results_gaps["Multi2B1"] = gap_m2b1;
            }
            {
                mdl.Solve(multicuts::Multi2B2);
                obj_m2b2.push_back(mdl.ObjectiveValue);
                time_m2b2.push_back(mdl.SolveTime);
                gap_m2b2.push_back(mdl.OptimalityGap);
                Results_times["Multi2B2"] = time_m2b2;
                Results_objs["Multi2B2"] = obj_m2b2;
                Results_gaps["Multi2B2"] = gap_m2b2;
            }
            {
                mdl.Solve(multicuts::Multi3B1);
                obj_m3b1.push_back(mdl.ObjectiveValue);
                time_m3b1.push_back(mdl.SolveTime);
                gap_m3b1.push_back(mdl.OptimalityGap);
                Results_times["Multi3B1"] = time_m3b1;
                Results_objs["Multi3B1"] = obj_m3b1;
                Results_gaps["Multi3B1"] = gap_m3b1;
            }
            {
                mdl.Solve(multicuts::Multi3B2);
                obj_m3b2.push_back(mdl.ObjectiveValue);
                time_m3b2.push_back(mdl.SolveTime);
                gap_m3b2.push_back(mdl.OptimalityGap);
                Results_times["Multi3B2"] = time_m3b2;
                Results_objs["Multi3B2"] = obj_m3b2;
                Results_gaps["Multi3B2"] = gap_m3b2;
            }
            cout << "\n" << endl;

            ///Important: reset the writing of the files once memory leak has been found and fixed
            std::ofstream f1(resultspath + "ObjectiveValues.json");
            f1 << std::setw(3) << Results_objs << std::endl;
            f1.close();
            std::ofstream f2(resultspath + "SolveTimes.json");
            f2 << std::setw(3) << Results_times << std::endl;
            f2.close();
            std::ofstream f3(resultspath + "OptimalityGaps.json");
            f3 << std::setw(3) << Results_gaps << std::endl;
            f3.close();
        }
        }
#endif

	return 0;
}

