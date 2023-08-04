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
     
    final["Solve time, MIP (sec)"] = (double)temp.value("Solve time, MIP (x100)", -1) / 100;
    final["Objective value"] = (double)temp.value("ObjectiveValue (x100)", -1) / 100;

    if (temp.contains("Solve time, LP (x100)")) {
        final["Solve time, LP (sec)"] = (double)temp.value("Solve time, LP (x100)", -1) / 100;
    }
    if (temp.contains("CplexStatus")) {
        final["Cplex status"] = StatusConversion[temp.value("CplexStatus", 0)];
    }
    if (temp.contains("OptimalityGap (x100)")) {
        final["Optimality gap (%)"] = (double)temp.value("OptimalityGap (x100)", -1) / 100;
    }
    if (temp.contains("nNodes")) {
        final["Number of nodes"] = temp.value("nNodes", -1);
    }
    if (temp.contains("nLazyCuts")) {
        final["Number of lazy cuts"] = temp.value("nLazyCuts", -1);
    }
    if (temp.contains("LazyCutTime (x1000)")) {
        final["Average lazy cut time (sec)"] = (double)temp.value("LazyCutTime (x1000)", -1) / 1000;
    }
    if (temp.contains("nUserCuts")) {
        final["Number of user cuts"] = temp.value("nUserCuts", -1);
    }
    if (temp.contains("UserCutTime (x1000)")) {
        final["Average user cut time (sec)"] = (double)temp.value("UserCutTime (x1000)", -1) / 1000;
    }
    if (temp.contains("nRestricted")) {
        final["Number of restricted subproblems"] = (int)temp.value("nRestricted", -1);
    }
    if (temp.contains("nDiversified")) {
        final["Number of diversified subproblems"] = (int)temp.value("nDiversified", -1);
    }
    if (temp.contains("nManualBranches")) {
        final["Number of local branching separations"] = (int)temp.value("nManualBranches", -1);
    }
    if (temp.contains("GRASP cut time (x100)")) {
        final["GRASP cut time (sec)"] = (double)temp.value("GRASP cut time (x100)", -1) / 100;
    }
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
    vector<string> doubles = { "Solve time (sec)" , "Objective value" , "Optimality gap (%)" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
    vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts" };
    vector<string> strings = { "Cplex status" };

    for (string key : doubles) { temp[key] = vector<double>(); }
    for (string key : ints) { temp[key] = vector<int>(); }
    for (string key : strings) { temp[key] = vector<string>(); }


    time_t start;
    time(&start);
    std::string basefile = "C:\\Users\\dobby\\Desktop\\Git Hub repo\\Charging-Station_Cpp\\";
    string sharedFile = basefile + "Shared_Home.json";
    string instanceFile = basefile + "MC0_Home_compressed.json";
    Data data;
    //data.load_fromUtilities(sharedFile, instanceFile, true, priceProfile::OnlyLevel3_Unperturbed);
    data.load(sharedFile, instanceFile, true);
    data.params["Stations_maxNewOutletsPerTimePeriod"] = 4;
    data.params["Stations_maxNewStationsPerTimePeriod"] = 2;
    //data.T = 2;

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



    {
        vector<float> results_cutting;
        vector<float> results_multicut;

        //Define constants for easier reading and writing
        int& T = data.T;
        int& M_bar = data.M_bar;
        vector<int>& P = data.P;

        time_t start;
        time(&start);
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);


        //Set parameters
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::TimeLimit, 7200);
        cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
        cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 2);
        cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
        cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
        cplex.setParam(IloCplex::Param::WorkMem, 100000);
        cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);

        cplex.setParam(IloCplex::Param::MIP::Display, 0);
        cplex.setParam(IloCplex::Param::Tune::Display, 0);
        cplex.setParam(IloCplex::Param::Simplex::Display, 0);
        cplex.setParam(IloCplex::Param::Sifting::Display, 0);
        cplex.setParam(IloCplex::Param::ParamDisplay, 0);
        cplex.setParam(IloCplex::Param::Network::Display, 0);
        cplex.setParam(IloCplex::Param::Conflict::Display, 0);


        //Create variables
        BoolVar2D x(env, T);
        for (int t = 0; t < T; t++) {
            x[t] = IloBoolVarArray(env, M_bar);
        }


        IloNumVarArray  theta(env, T, 0.0, IloInfinity);

        IloNumVar theta_obj(env);

        BoolVar2D trust(env, 256);
        for (int i = 0; i < 256; i++) {
            trust[i] = IloBoolVarArray(env, T);
        }


        //Constraints
        ////Budget, year 0
        IloExpr budget0(env);
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = data.params["station_coord"][j_bar][0];
            int k = data.params["station_coord"][j_bar][1];
            budget0 += (double)data.params["c"][0][j][k] * (x[0][j_bar] - (int)data.params["x0"][j][k]);
        }
        model.add(budget0 <= (double)data.params["B"][0]);
        budget0.end();

        ////Budget, year 1+
        for (int t = 1; t < T; t++) {
            IloExpr budget(env);
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                int j = data.params["station_coord"][j_bar][0];
                int k = data.params["station_coord"][j_bar][1];
                budget += (double)data.params["c"][t][j][k] * (x[t][j_bar] - x[t - 1][j_bar]);
            }
            model.add(budget <= (double)data.params["B"][t]);
            budget.end();
        }


        ////Can't remove outlets, year 0
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = data.params["station_coord"][j_bar][0];
            int k = data.params["station_coord"][j_bar][1];
            model.add(x[0][j_bar] >= (int)data.params["x0"][j][k]);
        }

        ////Can't remove outlets, year 1+
        for (int t = 1; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                model.add(x[t][j_bar] >= x[t - 1][j_bar]);
            }
        }

        ////At least k outlets, year 0+
        for (int t = 0; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                int k = data.params["station_coord"][j_bar][1];
                if (k > 1) {
                    model.add(x[t][j_bar] <= x[t][j_bar - 1]);
                }
            }
        }

        ///Redundant, but prevents crash from unextracted variables
        model.add(theta_obj >= 0);
        for (int t = 0; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                model.add(x[t][j_bar] >= 0);
            }
            model.add(theta[t] >= 0);
            for (auto i = 0; i < trust.getSize(); i++) {
                model.add(trust[i][t] >= 0);
            }
        }

        ////Upper bound for theta (ensures bounded problem)

        for (int t = 0; t < T; t++) {
            IloNum ub = data.weights[t].sum();
            model.add(theta[t] <= ub);
        }

        ////Ensures trust cut big-Ms work correctly
        ////Without these constraints, it's unclear how CPLEX handles the trust cuts
        for (auto i = 0; i < trust.getSize(); i++) {
            model.add(IloSum(trust[i]) <= T - 1);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        //Objective
        IloExpr obj(env);
        for (int t = 0; t < T; t++) {
            obj += theta[t] + data.Precovered[t];
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                obj += data.Ps[t](j_bar) * x[t][j_bar];
            }
        }
        //model.add(IloMaximize(env, obj));
        //obj.end();
        model.add(theta_obj == obj);
        obj.end();

        model.add(IloMaximize(env, theta_obj));

        random_selector<> choice{};
        for (int _i = 0; _i < 10; _i++) {

            //Create random, feasible solution
            ArrayXXd sol = ArrayXXd::Constant(T, data.params["M"], 0.0);
            for (int t = 0; t < T; t++) {
                double budget = data.params["B"][t];
                while (budget > 0) {
                    vector<int> candidates;
                    for (int j = 0; j < data.params["M"]; j++) {
                        if ((sol(t, j) < data.params["Mj"][j] - 1) && (data.params["c"][t][j][sol(t, j) + 1] <= budget)) { candidates.push_back(j); }
                    }
                    if (candidates.size() == 0) { break; }
                    int j_bar = choice(candidates);
                    if (j_bar == -1) { break;}
                    else {
                        for (int t_bar = t; t_bar < T; t_bar++) {
                            sol(t_bar, j_bar) += 1;
                        }
                        budget -= data.params["c"][t][j_bar][sol(t, j_bar)];
                    }

                }
            }
            sol = sol.round();
            cout << "Initial solution: ";
            for (int t = 0; t < T; t++) {
                cout << sol.row(t) << endl;
            }
            ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
            for (int t = 0; t < T; t++) {
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    int j = data.params["station_coord"][j_bar][0];
                    int k = data.params["station_coord"][j_bar][1];
                    if (sol(t, j) >= k) { x_tilde(t, j_bar) = 1; }
                }
            }
            vector<VectorXd> I_tilde;
            for (int t = 0; t < T; t++) {
                I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (x_tilde(t, j_bar) > 1 - EPS) {
                        I_tilde[t] += data.a[t].col(j_bar);
                    }
                }

            }
            double obj = 0.0;
            for (int t = 0; t < T; t++) {
                double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                obj += covered + data.Precovered[t];
                for (int j_bar = 0; j_bar < M_bar; j_bar++) { if (x_tilde(t, j_bar) > 1 - EPS) { obj += data.Ps[t](j_bar); } }
            }
            cout << "Objective for x_tilde: " << obj << endl;
            



            //Subproblem, multicut
            {
                cout << "Multicut method" << endl;
                chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
                IloModel submodel(env);
                submodel.add(model);

                IloCplex subcplex(submodel);
                subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 0);
                subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

                subcplex.setParam(IloCplex::Param::Threads, 1);
                subcplex.setParam(IloCplex::Param::Tune::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
                subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
                subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
                subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
                subcplex.setParam(IloCplex::Param::Network::Display, 0);
                subcplex.setParam(IloCplex::Param::Conflict::Display, 0);


                for (int t = 0; t < T; t++) {
                    IloExpr delta(env);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                        else { delta += x[t][j_bar]; }
                    }
                    submodel.add(delta <= 2.0);
                    delta.end();
                }
                IloNumVarArray startVar(env);
                IloNumArray startVal(env);
                for (int t = 0; t < T; t++) {
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        startVar.add(x[t][j_bar]);
                        startVal.add(x_tilde(t, j_bar));
                    }
                }
                subcplex.addMIPStart(startVar, startVal);
                startVal.end();
                startVar.end();

                ////Link callback
                MultiCutBenders_Callback cb(data, x, theta);
                CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
                    | IloCplex::Callback::Context::Id::Relaxation;
                subcplex.use(&cb, contextmask);


                bool solved = subcplex.solve();
                float time = (float)chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
                results_multicut.push_back(time);
                cout << "Time: " << time << endl;
                cout << "Objective: " << subcplex.getObjValue() << endl;
                ArrayXXd solution = ArrayXXd::Constant(T, data.params["M"], 0.0);
                for (int t = 0; t < T; t++) {
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        int j = data.params["station_coord"][j_bar][0];
                        int k = data.params["station_coord"][j_bar][1];
                        int val = subcplex.getValue(x[t][j_bar]);
                        if ((val > 1 - EPS) && (solution(t, j) < k)) { solution(t, j) = k; }
                    }
                }
                cout << "x_star \t theta \t obj " << endl;
                for (int t = 0; t < T; t++) {
                    cout << solution.row(t) << "\t";
                    cout << subcplex.getValue(theta[t]) << "\t";
                    VectorXd I_star = VectorXd::Constant(P[t], 0.0);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (subcplex.getValue(x[t][j_bar]) > 1 - EPS) {
                            I_star += data.a[t].col(j_bar);
                        }
                    }

                    double obj = 0.0;
                    double covered = (I_star.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                    obj += covered + data.Precovered[t];
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) { if (subcplex.getValue(x[t][j_bar]) > 1 - EPS) { obj += data.Ps[t](j_bar); } }
                    cout << covered << endl;
                }
                cout << endl;
                subcplex.end();
                submodel.end();
            }



            //Subproblem, cutting
            {
                cout << "Improved method" << endl;
                chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
                IloModel submodel(env);
                submodel.add(model);


                IloCplex subcplex(submodel);
                subcplex.use(NULL, 0);
                subcplex.setParam(IloCplex::Param::Threads, 1);
                subcplex.setParam(IloCplex::Param::Tune::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
                subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
                subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
                subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
                subcplex.setParam(IloCplex::Param::Network::Display, 0);
                subcplex.setParam(IloCplex::Param::Conflict::Display, 0);


                for (int t = 0; t < T; t++) {
                    IloExpr delta(env);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                        else { delta += x[t][j_bar]; }
                    }
                    submodel.add(delta <= 2.0);
                    delta.end();
                }

                //Two-opt cut generation
                for (int t = 0; t < T; t++) {
                    vector<int> active;
                    vector<int> inactive;
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { active.push_back(j_bar); }
                        else { inactive.push_back(j_bar); }

                    }
                    //Initialise coverage
                    VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
                    for (int j_bar : active) { I_tilde += data.a[t].col(j_bar); }

                    //Cut for greedy solution + distance 1
                    {
                        IloExpr lhs(env);
                        lhs -= theta[t];
                        double covered = (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
                        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                            lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        }
                        submodel.add(lhs >= -covered);
                    }

                    //Cuts for adding two outlets (redundant for greedy solution)
                    for (int j_bar : inactive) {
                        IloExpr lhs(env);
                        lhs -= theta[t];

                        VectorXd mod_I_tilde = I_tilde + data.a[t].col(j_bar);
                        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        for (int j_hat : inactive) {
                            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                        }
                        submodel.add(lhs >= -covered);
                    }

                    //Check coverage when swapping out active stations with different ones
                    for (int j_bar : active) {
                        IloExpr lhs(env);
                        lhs -= theta[t];

                        VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
                        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        for (int j_hat : inactive) {
                            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                        }
                        submodel.add(lhs >= -covered);
                    }
                }

                ////Two-opt cut generation
                //for (int t = 0; t < T; t++) {
                //    vector<int> always;
                //    vector<int> active;
                //    vector<int> inactiveD1;
                //    vector<int> inactiveD2;
                //    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                //        int k = data.params["station_coord"][j_bar][1];
                //        if (x_tilde(t, j_bar) > 1 - EPS) {
                //            int j = data.params["station_coord"][j_bar][0];
                //            if ((k == ((int)data.params["Mj"][j] - 1)) || (x_tilde(t, j_bar + 1) < EPS)) { active.push_back(j_bar); }
                //            else {
                //                always.push_back(j_bar);
                //                submodel.add(x[t][j_bar] == x_tilde(t, j_bar));
                //            }
                //        }
                //        else if ((k == 1) || (x_tilde(t, j_bar - 1) > 1 - EPS)) { inactiveD1.push_back(j_bar); }
                //        else if ((k == 2) || (x_tilde(t, j_bar - 2) > 1 - EPS)) { inactiveD2.push_back(j_bar); }
                //        else { submodel.add(x[t][j_bar] == x_tilde(t, j_bar)); }

                //    }
                //    //Initialise coverage
                //    VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
                //    for (int j_bar : always) { I_tilde += data.a[t].col(j_bar); }
                //    for (int j_bar : active) { I_tilde += data.a[t].col(j_bar); }

                //    //Cut for greedy solution + distance 1
                //    {
                //        IloExpr lhs(env);
                //        lhs -= theta[t];
                //        double covered = (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                //        VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
                //        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                //            lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                //        }
                //        submodel.add(lhs >= -covered);
                //    }

                //    //Cuts for adding two outlets (redundant for greedy solution)
                //    for (int j_bar : inactiveD1) {
                //        IloExpr lhs(env);
                //        lhs -= theta[t];

                //        VectorXd mod_I_tilde = I_tilde + data.a[t].col(j_bar);
                //        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                //        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                //        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                //        for (int j_hat : inactiveD2) {
                //            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                //        }
                //        submodel.add(lhs >= -covered);
                //    }

                //    //Check coverage when swapping out active stations with different ones
                //    for (int j_bar : active) {
                //        IloExpr lhs(env);
                //        lhs -= theta[t];

                //        VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
                //        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                //        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                //        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                //        for (int j_hat : inactiveD1) {
                //            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                //        }
                //        submodel.add(lhs >= -covered);
                //    }
                //}
                bool solved = subcplex.solve();
                float time = (float)chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
                results_cutting.push_back(time);
                cout << "Time: " << time << endl;
                cout << "Objective: " << subcplex.getObjValue() << endl;             
                ArrayXXd solution = ArrayXXd::Constant(T, data.params["M"], 0.0);
                for (int t = 0; t < T; t++) {
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        int j = data.params["station_coord"][j_bar][0];
                        int k = data.params["station_coord"][j_bar][1];
                        int val = subcplex.getValue(x[t][j_bar]);
                        if ((val > 1 - EPS) && (solution(t,j) < k )){ solution(t, j) = k; }
                    }
                }
                cout << "x_star \t theta \t obj " << endl;
                for (int t = 0; t < T; t++) {
                    cout << solution.row(t) << "\t";
                    cout << subcplex.getValue(theta[t]) << "\t";
                    VectorXd I_star = VectorXd::Constant(P[t], 0.0);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (subcplex.getValue(x[t][j_bar]) > 1 - EPS) {
                            I_star += data.a[t].col(j_bar);
                        }
                    }

                    double obj = 0.0;
                    double covered = (I_star.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                    obj += covered + data.Precovered[t];
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) { if (subcplex.getValue(x[t][j_bar]) > 1 - EPS) { obj += data.Ps[t](j_bar); } }
                    cout << covered << endl;
                }

                cout << endl;
                subcplex.end();
                submodel.end();
            }


            ////Add Benders optimality cuts and local branching cuts
            //for (int t = 0; t < T; t++) {
            //    IloExpr lhs(env);
            //    lhs -= theta[t];

            //    //Multi-cut (by time period), B1
            //    double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
            //    VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
            //    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            //        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
            //    }
            //    model.add(lhs >= -covered);
            //}

            for (int t = 0; t < T; t++) {
                IloExpr delta(env);
                delta += 3.0 * trust[_i][t];
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                    else { delta += x[t][j_bar]; }
                }
                model.add(delta >= 3.0);
            }

        }





        env.end();
    }





#endif

#ifdef __linux__
    int maxTest = 20;
    std::string file;
    string dataset = argv[1];
    {
        cout << "Dataset: " << dataset << endl;
        std::string resultspath = "/local_1/outer/lamste/Results/TroisRivieres/C++/" + dataset + "/";
        string prefix = "SubproblemComparison";
        std::string fp_stats = resultspath + prefix + "_Statistics.json";

        
        std::string basefile = "/local_1/outer/lamste/Data/Precomputed/" + dataset + "/MaximumCover/";
        string sharedFile = basefile + "Shared.json";
        json Results_stats = { {"Multicut", json()}, {"Cutting", json()} };


        Results_stats["Multicut"]["Times"] = vector<vector<float>>();
        Results_stats["Multicut"]["Objectives"] = vector<vector<float>>();
        Results_stats["Cutting"]["Times"] = vector<vector<float>>();
        Results_stats["Cutting"]["Objectives"] = vector<vector<float>>();
        //{
        //    json Results;
        //    std::ifstream f(fp_stats, ifstream::in);
        //    f >> Results;
        //    f.close();
        //    Results_stats = Results;
        //}
        //vector<string> doubles = { "Solve time, LP (sec)" , "Solve time, MIP (sec)", "Objective value" , "Optimality gap (%)" , "Average lazy cut time (sec)" , "Average user cut time (sec)" };
        //vector<string> ints = { "Number of nodes" , "Number of lazy cuts" , "Number of user cuts", "Number of restricted subproblems", "Number of diversified subproblems", "Number of local branching separations" };
        //vector<string> strings = { "Cplex status" };
        //vector<string> labels = {"LocalBranching_MoreBranching"};
        //for (string label : labels) {
        //    if (Results_stats.contains(label)){ continue; }
        //    Results_stats[label] = {};
        //    for (string key : doubles) { Results_stats[label][key] = vector<double>(); }
        //    for (string key : ints) { Results_stats[label][key] = vector<int>(); }
        //    for (string key : strings) { Results_stats[label][key] = vector<string>(); }
        //}
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

            //data.T = 2; //Remember to remove this later
            //data.load_fromUtilities(sharedFile, instanceFile, true, priceProfile::OnlyLevel3_Unperturbed);

            cout << "Data loading time: " << time(NULL) - start << " seconds" << endl;
            cout << endl;
            
                        
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

            //{
            //    string label = "LocalBranching_MoreBranching";
            //    if (!(Results_stats.contains(label))) { throw std::runtime_error("Label missing from JSON. Add label to labels vector."); }
            //    if (Results_stats[label]["Cplex status"].size() > (long unsigned int) test) { throw std::runtime_error("Results already populated. Check the label is spelled correctly or disable this error."); }
            //    json params = { { "verbose", true }, {"budgetType", BUDGET_TYPE::Knapsack} };
            //    LocalBranching_Model mdl;
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

        {
        vector<float> times_cutting;
        vector<float> objs_cutting;
        vector<float> times_multicut;
        vector<float> objs_multicut;


        //Define constants for easier reading and writing
        int& T = data.T;
        int& M_bar = data.M_bar;
        vector<int>& P = data.P;

        time_t start;
        time(&start);
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);


        //Set parameters
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::TimeLimit, 7200);
        cplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 2);
        cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 2);
        cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
        cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
        cplex.setParam(IloCplex::Param::WorkMem, 100000);
        cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);

        cplex.setParam(IloCplex::Param::MIP::Display, 0);
        cplex.setParam(IloCplex::Param::Tune::Display, 0);
        cplex.setParam(IloCplex::Param::Simplex::Display, 0);
        cplex.setParam(IloCplex::Param::Sifting::Display, 0);
        cplex.setParam(IloCplex::Param::ParamDisplay, 0);
        cplex.setParam(IloCplex::Param::Network::Display, 0);
        cplex.setParam(IloCplex::Param::Conflict::Display, 0);


        //Create variables
        BoolVar2D x(env, T);
        for (int t = 0; t < T; t++) {
            x[t] = IloBoolVarArray(env, M_bar);
        }


        IloNumVarArray  theta(env, T, 0.0, IloInfinity);

        IloNumVar theta_obj(env);

        BoolVar2D trust(env, 256);
        for (int i = 0; i < 256; i++) {
            trust[i] = IloBoolVarArray(env, T);
        }


        //Constraints
        ////Budget, year 0
        IloExpr budget0(env);
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = data.params["station_coord"][j_bar][0];
            int k = data.params["station_coord"][j_bar][1];
            budget0 += (double)data.params["c"][0][j][k] * (x[0][j_bar] - (int)data.params["x0"][j][k]);
        }
        model.add(budget0 <= (double)data.params["B"][0]);
        budget0.end();

        ////Budget, year 1+
        for (int t = 1; t < T; t++) {
            IloExpr budget(env);
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                int j = data.params["station_coord"][j_bar][0];
                int k = data.params["station_coord"][j_bar][1];
                budget += (double)data.params["c"][t][j][k] * (x[t][j_bar] - x[t - 1][j_bar]);
            }
            model.add(budget <= (double)data.params["B"][t]);
            budget.end();
        }


        ////Can't remove outlets, year 0
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = data.params["station_coord"][j_bar][0];
            int k = data.params["station_coord"][j_bar][1];
            model.add(x[0][j_bar] >= (int)data.params["x0"][j][k]);
        }

        ////Can't remove outlets, year 1+
        for (int t = 1; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                model.add(x[t][j_bar] >= x[t - 1][j_bar]);
            }
        }

        ////At least k outlets, year 0+
        for (int t = 0; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                int k = data.params["station_coord"][j_bar][1];
                if (k > 1) {
                    model.add(x[t][j_bar] <= x[t][j_bar - 1]);
                }
            }
        }

        ///Redundant, but prevents crash from unextracted variables
        model.add(theta_obj >= 0);
        for (int t = 0; t < T; t++) {
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                model.add(x[t][j_bar] >= 0);
            }
            model.add(theta[t] >= 0);
        }

        ////Upper bound for theta (ensures bounded problem)

        for (int t = 0; t < T; t++) {
            IloNum ub = data.weights[t].sum() + data.Precovered[t] + data.Ps[t].sum();
            model.add(theta[t] <= ub);
        }

        ////Ensures trust cut big-Ms work correctly
        ////Without these constraints, it's unclear how CPLEX handles the trust cuts
        for (auto i = 0; i < trust.getSize(); i++) {
            model.add(IloSum(trust[i]) <= T - 1);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////

        //Objective
        IloExpr obj(env);
        for (int t = 0; t < T; t++) {
            obj += theta[t] + data.Precovered[t];
            for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                obj += data.Ps[t](j_bar) * x[t][j_bar];
            }
        }
        model.add(IloMaximize(env, obj));
        obj.end();

        random_selector<> choice{};
        for (int _i = 0; _i < 255; _i++) {

            //Create random, feasible solution
            ArrayXXd sol = ArrayXXd::Constant(T, data.params["M"], 0.0);
            for (int t = 0; t < T; t++) {
                double budget = data.params["B"][t];
                while (budget > 0) {
                    vector<int> candidates;
                    for (int j = 0; j < data.params["M"]; j++) {
                        if ((sol(t, j) < (int) data.params["Mj"][j] - 1) && ((double) data.params["c"][t][j][sol(t, j) + 1] <= budget)) { candidates.push_back(j); }
                    }
                    if (candidates.size() == 0) { break; }
                    int j_bar = choice(candidates);
                    if (j_bar == -1) { break; }
                    else {
                        for (int t_bar = t; t_bar < T; t_bar++) {
                            sol(t_bar, j_bar) += 1;
                        }
                        budget -= (double) data.params["c"][t][j_bar][sol(t, j_bar)];
                    }

                }
            }
            sol = sol.round();
            ArrayXXd x_tilde = ArrayXXd::Constant(T, M_bar, 0.0);
            for (int t = 0; t < T; t++) {
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    int j = data.params["station_coord"][j_bar][0];
                    int k = data.params["station_coord"][j_bar][1];
                    if (sol(t, j) >= k) { x_tilde(t, j_bar) = 1; }
                }
            }
            vector<VectorXd> I_tilde;
            for (int t = 0; t < T; t++) {
                I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (x_tilde(t, j_bar) > 1 - EPS) {
                        I_tilde[t] += data.a[t].col(j_bar);
                    }
                }

            }


            //Subproblem, multicut
            {
                chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
                IloModel submodel(env);
                submodel.add(model);

                IloCplex subcplex(submodel);
                subcplex.setParam(IloCplex::Param::Preprocessing::Reformulations, 0);
                subcplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

                subcplex.setParam(IloCplex::Param::Threads, 1);
                subcplex.setParam(IloCplex::Param::Tune::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
                subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
                subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
                subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
                subcplex.setParam(IloCplex::Param::Network::Display, 0);
                subcplex.setParam(IloCplex::Param::Conflict::Display, 0);
                subcplex.setParam(IloCplex::Param::TimeLimit, 60);


                for (int t = 0; t < T; t++) {
                    IloExpr delta(env);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                        else { delta += x[t][j_bar]; }
                    }
                    submodel.add(delta <= 2.0);
                    delta.end();
                }

                ////Link callback
                MultiCutBenders_Callback cb(data, x, theta);
                CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
                    | IloCplex::Callback::Context::Id::Relaxation;
                subcplex.use(&cb, contextmask);


                bool solved = subcplex.solve();
                float time = (float)chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
                times_multicut.push_back(time);
                objs_multicut.push_back(subcplex.getObjValue());
                subcplex.end();
                submodel.end();
            }

            //Subproblem, cutting
            {
                chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
                IloModel submodel(env);
                submodel.add(model);


                IloCplex subcplex(submodel);
                subcplex.setParam(IloCplex::Param::Threads, 1);
                subcplex.setParam(IloCplex::Param::Tune::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Display, 0);
                subcplex.setParam(IloCplex::Param::MIP::Interval, 10000);
                subcplex.setParam(IloCplex::Param::Simplex::Display, 0);
                subcplex.setParam(IloCplex::Param::Sifting::Display, 0);
                subcplex.setParam(IloCplex::Param::ParamDisplay, 0);
                subcplex.setParam(IloCplex::Param::Network::Display, 0);
                subcplex.setParam(IloCplex::Param::Conflict::Display, 0);
                subcplex.setParam(IloCplex::Param::TimeLimit, 60);


                for (int t = 0; t < T; t++) {
                    IloExpr delta(env);
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                        else { delta += x[t][j_bar]; }
                    }
                    submodel.add(delta <= 2.0);
                    delta.end();
                }

                //Two-opt cut generation
                for (int t = 0; t < T; t++) {
                    vector<int> active;
                    vector<int> inactive;
                    for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                        if (x_tilde(t, j_bar) > 1 - EPS) { active.push_back(j_bar); }
                        else { inactive.push_back(j_bar); }

                    }
                    //Initialise coverage
                    VectorXd I_tilde = VectorXd::Constant(P[t], 0.0);
                    for (int j_bar : active) { I_tilde += data.a[t].col(j_bar); }

                    //Cut for greedy solution + distance 1
                    {
                        IloExpr lhs(env);
                        lhs -= theta[t];
                        double covered = (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (I_tilde.array() < 1).matrix().cast<double>();
                        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                            lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        }
                        submodel.add(lhs >= -covered);
                    }

                    //Cuts for adding two outlets (redundant for greedy solution)
                    for (int j_bar : inactive) {
                        IloExpr lhs(env);
                        lhs -= theta[t];

                        VectorXd mod_I_tilde = I_tilde + data.a[t].col(j_bar);
                        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        for (int j_hat : inactive) {
                            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                        }
                        submodel.add(lhs >= -covered);
                    }

                    //Check coverage when swapping out active stations with different ones
                    for (int j_bar : active) {
                        IloExpr lhs(env);
                        lhs -= theta[t];

                        VectorXd mod_I_tilde = I_tilde - data.a[t].col(j_bar);
                        double covered = (mod_I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                        VectorXd uncovered = (mod_I_tilde.array() < 1).matrix().cast<double>();
                        lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                        for (int j_hat : inactive) {
                            lhs += (data.CutCoeffs[t].col(j_hat).dot(uncovered)) * x[t][j_hat];
                        }
                        submodel.add(lhs >= -covered);
                    }
                }
                bool solved = subcplex.solve();
                float time = (float)chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
                times_cutting.push_back(time);
                objs_cutting.push_back(subcplex.getObjValue());
                subcplex.end();
                submodel.end();
            }


            //Add Benders optimality cuts and local branching cuts
            for (int t = 0; t < T; t++) {
                IloExpr lhs(env);
                lhs -= theta[t];

                //Multi-cut (by time period), B1
                double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
                VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
                }
                model.add(lhs >= -covered);
            }

            for (int t = 0; t < T; t++) {
                IloExpr delta(env);
                delta += 3.0 * trust[_i][t];
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (x_tilde(t, j_bar) > 1 - EPS) { delta += 1 - x[t][j_bar]; }
                    else { delta += x[t][j_bar]; }
                }
                model.add(delta >= 3.0);
            }

        }
        Results_stats["Multicut"]["Times"].push_back(times_multicut);
        Results_stats["Multicut"]["Objectives"].push_back(objs_multicut);
        Results_stats["Cutting"]["Times"].push_back(times_cutting);
        Results_stats["Cutting"]["Objectives"].push_back(objs_cutting);

        env.end();
        }
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



