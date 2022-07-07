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

#include "Model_BendersCordeau.h"
#include "CoverageCallback.h"


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
    
	std::string file = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC0.json";
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
    vector<float> Results;
    vector<float> Times;

#ifdef __linux__
    std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/Simple/";
    std::string basefile = "/local_1/outer/lamste/Data/Precomputed/Simple/MaximumCover/MC%d.json";
    int maxTest = 20;
#endif
#ifdef _WIN32
    std::string resultspath = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\";
    std::string basefile =  "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\MC%d_Home.json";
    int maxTest = 1;
#endif

    for (int test = 0; test < maxTest; test++) {
        time_t start;
        time(&start);
        file = string_format(basefile, test);
        data.load(file, verbose);
        cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
        
        //vector<vector<int>> oldSol = {
        //    {
        //        0,
        //        0,
        //        0,
        //        2,
        //        2,
        //        0,
        //        0,
        //        0,
        //        0,
        //        0
        //    },
        //        {
        //            2,
        //            0,
        //            2,
        //            2,
        //            2,
        //            0,
        //            0,
        //            0,
        //            0,
        //            0
        //        },
        //        {
        //            2,
        //            2,
        //            2,
        //            2,
        //            2,
        //            0,
        //            0,
        //            0,
        //            2,
        //            0
        //        },
        //        {
        //            2,
        //            2,
        //            2,
        //            2,
        //            2,
        //            2,
        //            0,
        //            0,
        //            2,
        //            2
        //        }
        //};

        Model_BendersCordeau mdl;
        //Model_BB mdl;
        mdl.SetData(data);
        mdl.Solve();
        Results.push_back(mdl.ObjectiveValue);
        Times.push_back(mdl.SolveTime);


        json results;
        results["ObjectiveValue"] = Results;
        results["SolveTime"] = Times;
        results["Solution"] = mdl.Solution;
        std::ofstream f(resultspath + "results.json");
        f << std::setw(3) << results << std::endl;



    }

	return 0;
}

