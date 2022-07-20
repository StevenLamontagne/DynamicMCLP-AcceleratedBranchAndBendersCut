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

    //Greedy mdl;
    //mdl.SetData(data);
    //mdl.Solve();
    //cout << "Solution value (Greedy): " << mdl.SolutionQuality << endl;
    //cout << "Solution value (Data): " << data.SolutionQuality(mdl.Solution) << endl;
    //cout << "Solve time (seconds): " << mdl.SolveTime << endl;

    Model_Multicut mdl2;
    mdl2.SetData(data);
    mdl2.Solve(multicuts::Multi1B1);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;

    std::cout << "\n" << endl;
    std::cout << "\n" << endl;

    mdl2.Solve(multicuts::Multi1B2);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;

    std::cout << "\n" << endl;
    std::cout << "\n" << endl;

    mdl2.Solve(multicuts::Multi2B1);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;

    std::cout << "\n" << endl;
    std::cout << "\n" << endl;

    mdl2.Solve(multicuts::Multi2B2);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;

    std::cout << "\n" << endl;
    std::cout << "\n" << endl;

    mdl2.Solve(multicuts::Multi3B1);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;

    std::cout << "\n" << endl;
    std::cout << "\n" << endl;

    mdl2.Solve(multicuts::Multi3B2);
    std::cout << "Solution value (Benders): " << mdl2.ObjectiveValue << endl;
    std::cout << "Solution value (Data): " << data.SolutionQuality(mdl2.Solution) << endl;
    std::cout << "Solve time (seconds): " << mdl2.SolveTime << endl;
    //Results_times["CutTimes"] = *mdl2.CutTimes;
    //std::ofstream f2(resultspath + "CutTimes.json");
    //f2 << std::setw(3) << Results_times << std::endl;
    //f2.close();




#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC0.json";
    vector<std::string> datasets = {"Price", "LongSpan", "HomeCharging", };
    //vector<std::string> datasets = { "Distance" };


    for (int dBar = 0; dBar < datasets.size(); dBar++) {
        cout << "Dataset: " << datasets[dBar] << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + datasets[dBar] + "/";
        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + datasets[dBar] + "/MaximumCover/MC%d.json";
        json Results_times;
        json Results_objs;
        json Results_gaps;

        vector<double> obj_bb;
        vector<double> time_bb;
        vector<double> gap_bb;


        vector<double> obj_sb1;
        vector<double> time_sb1;
        vector<double> gap_sb1;


        vector<double> obj_sb2;
        vector<double> time_sb2;
        vector<double> gap_sb2;

        for (int test = 0; test < maxTest; test++) {
            cout << "Test: " << test << endl;
            time_t start;
            time(&start);
            file = string_format(basefile, test);
            data.load(file, verbose);
            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;

            {
                Model_BendersCordeau mdl_sb1;
                mdl_sb1.SetData(data);
                mdl_sb1.Solve(cuts::SingleB1);
                obj_sb1.push_back(mdl_sb1.ObjectiveValue);
                time_sb1.push_back(mdl_sb1.SolveTime);
                gap_sb1.push_back(mdl_sb1.OptimalityGap);
                Results_times["SingleB1"] = time_sb1;
                Results_objs["SingleB1"] = obj_sb1;
                Results_gaps["SingleB1"] = gap_sb1;
            }
            {
                Model_BendersCordeau mdl_sb2;
                mdl_sb2.SetData(data);
                mdl_sb2.Solve(cuts::SingleB2);
                obj_sb2.push_back(mdl_sb2.ObjectiveValue);
                time_sb2.push_back(mdl_sb2.SolveTime);
                gap_sb2.push_back(mdl_sb2.OptimalityGap);
                Results_times["SingleB2"] = time_sb2;
                Results_objs["SingleB2"] = obj_sb2;
                Results_gaps["SingleB2"] = gap_sb2;
            }

            {
                Model_BB mdl_bb;
                mdl_bb.SetData(data);
                mdl_bb.Solve();
                obj_bb.push_back(mdl_bb.ObjectiveValue);
                time_bb.push_back(mdl_bb.SolveTime);
                gap_bb.push_back(mdl_bb.OptimalityGap);
                Results_times["B&B"] = time_bb;
                Results_objs["B&B"] = obj_bb;
                Results_gaps["B&B"] = gap_bb;
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

