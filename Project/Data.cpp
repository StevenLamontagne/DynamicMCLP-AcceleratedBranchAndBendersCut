#include "Data.h"



void Data::load(string sharedFile, string instanceFile, bool verbose)
{
    if (verbose) {
        cout << "File loading starting \n";
    }
    std::ifstream f(sharedFile, ifstream::in);

    try
    {
        f >> params;
        T = (int)params["T"];
        

        if (verbose) {
            cout << "File loading success \n";
            cout << "Covering creation started \n";
        }
        try {
            create_covering(instanceFile, verbose);
        }
        catch (exception & e) {
            cout << "Exception: " << e.what() << endl;
        }
        if (verbose) {
            cout << "Covering creation successful \n";
        }
    }
    catch (json::exception & e)
    {
        cout << e.what() << '\n';
    }

    params["Stations_isLevel3"] = isLevel3;
    params["Stations_costMultiplier"] = costMultiplier;

    params["Stations_maxNewOutletsPerTimePeriod"] = 4;
    params["Stations_maxNewStationsPerTimePeriod"] = 2;


    f.close();


}

void Data::load_fromUtilities(string sharedFile, string instanceFile, bool verbose, priceProfile priceProfile)
{
    if (verbose) {
        cout << "File loading starting \n";
    }
    std::ifstream f(sharedFile, ifstream::in);


    f >> params;
    T = (int)params["T"];

    a.clear();
    Ps.clear();
    Precovered.clear();
    M_bar = (int)params["station_coord"].size();
    int P_est = (int)params["user_coord"].size();
    for (int t = 0; t < T; t++) {
        Ps.push_back(Eigen::VectorXd::Zero(M_bar));
        Precovered.push_back(0.0);
    }

    if (verbose) {
        cout << "Loading coverage file \n";
    }

    params["Stations_isLevel3"] = isLevel3;
    params["Stations_costMultiplier"] = costMultiplier;


    int maxOutlets = 0;
    for (int j = 0; j < params["M"]; j++) {
        maxOutlets += (int) params["Mj"][j] - 1;
    }
    params["Stations_maxNewOutletsPerTimePeriod"] = 4;
    params["Stations_maxNewStationsPerTimePeriod"] = 2;


    json instance;
    std::ifstream f2(instanceFile, ifstream::in);
    f2 >> instance;
    f2.close();


    
    vector<vector<vector<double>>> c;
    for (int t = 0; t < T; t++) {
        vector<vector<double>> c1;
        for (int j = 0; j < params["M"]; j++) {
            vector<double> c2;
            c2.push_back(0.0);
            
            for (int k = 1; k < params["Mj"][j]; k++) {

                //Use different price profiles for level 2 versus level 3
                switch (priceProfile)
                {
                case priceProfile::MixedLevel2andLevel3_Unperturbed:
                {
                    if (isLevel3[j]) {
                        double val = 0.0;
                        if (k % 4 == 1) { val = (double)instance["c"][t][j] + (double)instance["f"][t][j]; }
                        else { val = (double)instance["c"][t][j]; }
                        c2.push_back(val);
                    }
                    else {
                        double val = 0.0;
                        if (k == 1) { val = 0.1 * (double)instance["c"][t][j] + 0.25 * (double)instance["f"][t][j]; }
                        else { val = 0.1 * (double)instance["c"][t][j]; }
                        c2.push_back(val);
                    }
                    break;
                }
                case priceProfile::MixedLevel2andLevel3_Perturbed:
                {
                    if (isLevel3[j]) {
                        double val = 0.0;
                        if (k % 4 == 1) { val = (double)instance["c"][t][j] + (double)instance["f"][t][j]; }
                        else { val = (double)instance["c"][t][j]; }
                        val *= costMultiplier[j]; 
                        c2.push_back(val);
                    }
                    else {
                        double val = 0.0;
                        if (k == 1) { val = 0.1 * (double)instance["c"][t][j] + 0.25 * (double)instance["f"][t][j]; }
                        else { val = 0.1 * (double)instance["c"][t][j]; }
                        val *= costMultiplier[j]; 
                        c2.push_back(val);
                    }
                    break;
                }
                case priceProfile::OnlyLevel3_Perturbed:
                {
                    double val = 0.0;
                    if (k == 1) { val = (double)instance["c"][t][j] + (double)instance["f"][t][j]; }
                    else { val = (double)instance["c"][t][j]; }
                    val *= costMultiplier[j];
                    c2.push_back(val);
                    break;
                }
                case priceProfile::OnlyLevel3_Unperturbed:
                default:
                {
                    double val = 0.0;
                    if (k == 1) { val = (double)instance["c"][t][j] + (double)instance["f"][t][j]; }
                    else { val = (double)instance["c"][t][j]; }
                    c2.push_back(val);
                    break;
                }
                }

            }
            c1.push_back(c2);
        }
        c.push_back(c1);
    }
    params["c"] = c;

    vector<vector<bool>> x0;
    for (int j = 0; j < params["M"]; j++) {
        vector<bool> x0_1;
        x0_1.push_back(0);
        for (int k = 1; k < params["Mj"][j]; k++) {
            if (instance["x0"][j] >= k) { x0_1.push_back(1); }
            else { x0_1.push_back(0); }
        }
        x0.push_back(x0_1);
    }
    params["x0"] = x0;


    //Calculate home coverage and station coverage
    for (int t = 0; t < T; t++) {

        Eigen::Array<bool, Eigen::Dynamic, 1> temp_home = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(P_est, 0);
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> temp_a = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(P_est, M_bar, 0);

        vector<double> temp_weight;
        std::vector<Tripletd> tripletList_a;
        std::vector<Tripletd> tripletList_CutCoeffs;
        tripletList_a.reserve(P_est);
        tripletList_CutCoeffs.reserve(P_est);


        for (auto p = 0; p < params["user_coord"].size(); p++) {
            int i = params["user_coord"][p][0];
            int r = params["user_coord"][p][1];
            double optout = instance["d0"][t][0][i][r];
            if ((instance["d0"][t][1][i][r] != nullptr) && (instance["d0"][t][1][i][r] >= optout) ){
                temp_home(p) = 1;
            }

            for (auto j_bar = 0; j_bar < params["station_coord"].size(); j_bar++) {
                int j = params["station_coord"][j_bar][0];
                int k = params["station_coord"][j_bar][1];
                if ((k > 1) && (temp_a(p, j_bar - 1))) {continue;}
                if (instance["d1"][t][j][i][r] != nullptr) {
                    double d1 = instance["d1"][t][j][i][r]; 
                    if (((priceProfile == priceProfile::MixedLevel2andLevel3_Perturbed) || (priceProfile == priceProfile::MixedLevel2andLevel3_Unperturbed)) && (!isLevel3[j]))
                    { d1 -= 1.463601; } //Reduce the value here if the station is level 2 instead of level 3 (by 1.463601, the value calculated by Mahsa)
                    if (instance["beta"][t][j][i][k] == nullptr) { continue; }
                    double beta = instance["beta"][t][j][i][k];
                    double u = beta + d1;
                    if (u >= optout) { temp_a(p, j_bar) = 1; }
                }
            }

            double weight = (double)params["Ni"][t][i] / (double)params["R"][i];

            if (temp_home(p) == 1) { Precovered[t] += weight; continue; }

            Eigen::Array<bool, 1, Eigen::Dynamic> sub = temp_a.row(p);
            int n = sub.count();
            switch (n)
            {
            case 0:
                break;
            case 1:
            {
                int j_bar;
                sub.maxCoeff(&j_bar);
                int j = params["station_coord"][j_bar][0];
                int k = params["station_coord"][j_bar][1];
                if (params["x0"][j][k] == 1) { Precovered[t] += weight; }
                else { Ps[t](j_bar) += weight; }
                break;
            }
            default:
            {
                bool pre = false;
                int reindex = temp_weight.size();
                vector<Tripletd> temp_Tripleta;
                vector<Tripletd> temp_TripletCC;
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (temp_a(p, j_bar)) {
                        int j = params["station_coord"][j_bar][0];
                        int k = params["station_coord"][j_bar][1];
                        if (params["x0"][j][k] == 1) { Precovered[t] += weight; pre = true; break; }
                        else { temp_Tripleta.push_back(Tripletd(reindex, j_bar, 1.0)); temp_TripletCC.push_back(Tripletd(reindex, j_bar, weight)); }
                    }
                }
                if (!pre) {
                    tripletList_a.insert(tripletList_a.end(), temp_Tripleta.begin(), temp_Tripleta.end());
                    tripletList_CutCoeffs.insert(tripletList_CutCoeffs.end(), temp_TripletCC.begin(), temp_TripletCC.end());
                    temp_weight.push_back(weight);
                }
                break;
            }
            }

        }

        P.push_back((int)temp_weight.size());
        Eigen::VectorXd convert = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(temp_weight.data(), temp_weight.size());
        weights.push_back(convert);

        SparseXXd new_a(P[t], M_bar);
        new_a.setFromTriplets(tripletList_a.begin(), tripletList_a.end());
        new_a.makeCompressed();
        a.push_back(new_a);


        SparseXXd new_CutCoeffs(P[t], M_bar);
        new_CutCoeffs.setFromTriplets(tripletList_CutCoeffs.begin(), tripletList_CutCoeffs.end());
        new_CutCoeffs.makeCompressed();
        CutCoeffs.push_back(new_CutCoeffs);

    }




    f.close();
}

double Data::SolutionQuality(vector<vector<int>> Sol)
{
    double total = 0;
    for (int t = 0; t < T; t++) {
        VectorXd cover = VectorXd::Constant(P[t], 0.0);
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = params["station_coord"][j_bar][0];
            int k = params["station_coord"][j_bar][1];
            if (Sol[t][j] >= k) { 
                cover += a[t].col(j_bar); 
                total += Ps[t](j_bar); }
        }
        total += weights[t].dot((VectorXd)(cover.array() >= 1).matrix().cast<double>());
        total += Precovered[t];
    }
    return total;
}

double Data::SolutionQuality(ArrayXXd Sol)
{
    double total = 0;
    for (int t = 0; t < T; t++) {
        VectorXd cover = VectorXd::Constant(P[t], 0.0);
        for (int j_bar = 0; j_bar < M_bar; j_bar++) {
            int j = params["station_coord"][j_bar][0];
            int k = params["station_coord"][j_bar][1];
            if (Sol(t,j) >= k) {
                cover += a[t].col(j_bar);
                total += Ps[t](j_bar);
            }
        }
        total += weights[t].dot((VectorXd)(cover.array() >= 1).matrix().cast<double>());
        total += Precovered[t];
    }
    return total;
}


void Data::create_covering(string instanceFile, bool verbose)
{
    a.clear();
    Ps.clear();
    Precovered.clear();

    if (verbose) {
        cout << "Loading coverage file \n";
    }


    json temp_coverage;
    std::ifstream f2(instanceFile, ifstream::in);
    f2 >> temp_coverage;
    f2.close();

    M_bar = (int) params["station_coord"].size();
    int P_est = (int) params["user_coord"].size();
    for (int t = 0; t < T; t++) {
        Ps.push_back(Eigen::VectorXd::Zero(M_bar));
        Precovered.push_back(0.0);
    }


    for (int t = 0; t < T; t++) {
        Eigen::Array<bool, Eigen::Dynamic, 1> temp_home = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(P_est, 0);
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> temp_a = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(P_est, M_bar, 0);

        vector<double> temp_weight;
        std::vector<Tripletd> tripletList_a;
        std::vector<Tripletd> tripletList_CutCoeffs;
        tripletList_a.reserve(P_est);
        tripletList_CutCoeffs.reserve(P_est);

        for (int coord : temp_coverage["Home"][t]) {
            temp_home(coord) = 1;
        }

        for (vector<int> coord : temp_coverage["a"][t]) {
            temp_a(coord[0],coord[1]) = 1;
        }


        for (int coord = 0; coord < P_est; coord++) {

            int i = params["user_coord"][coord][0];
            int r = params["user_coord"][coord][1];
            double weight = (double)params["Ni"][t][i] / (double)params["R"][i];

            if (temp_home(coord) == 1) { Precovered[t] += weight; continue; }

            Eigen::Array<bool, 1, Eigen::Dynamic> sub = temp_a.row(coord);
            int n = sub.count();
            switch (n)
            {
            case 0:
                break;
            case 1:
            {
                int j_bar;
                sub.maxCoeff(&j_bar);
                int j = params["station_coord"][j_bar][0];
                int k = params["station_coord"][j_bar][1];
                if (params["x0"][j][k] == 1) { Precovered[t] += weight; }
                else { Ps[t](j_bar) += weight; }
                break;
            }
            default:
            {
                bool pre = false;
                int reindex = temp_weight.size();
                vector<Tripletd> temp_Tripleta;
                vector<Tripletd> temp_TripletCC;
                for (int j_bar = 0; j_bar < M_bar; j_bar++) {
                    if (temp_a(coord, j_bar)) {
                        int j = params["station_coord"][j_bar][0];
                        int k = params["station_coord"][j_bar][1];
                        if (params["x0"][j][k] == 1) { Precovered[t] += weight; pre = true; break; }
                        else { temp_Tripleta.push_back(Tripletd(reindex, j_bar, 1.0)); temp_TripletCC.push_back(Tripletd(reindex, j_bar, weight)); }
                    }
                }
                if (!pre) {
                    tripletList_a.insert(tripletList_a.end(), temp_Tripleta.begin(), temp_Tripleta.end());
                    tripletList_CutCoeffs.insert(tripletList_CutCoeffs.end(), temp_TripletCC.begin(), temp_TripletCC.end());
                    temp_weight.push_back(weight);
                }
                break;
            }
            }
        }
        
        P.push_back((int) temp_weight.size());
        Eigen::VectorXd convert = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(temp_weight.data(), temp_weight.size());
        weights.push_back(convert);

        SparseXXd new_a(P[t], M_bar);
        new_a.setFromTriplets(tripletList_a.begin(), tripletList_a.end());
        new_a.makeCompressed();
        a.push_back(new_a);


        SparseXXd new_CutCoeffs(P[t], M_bar);
        new_CutCoeffs.setFromTriplets(tripletList_CutCoeffs.begin(), tripletList_CutCoeffs.end());
        new_CutCoeffs.makeCompressed();
        CutCoeffs.push_back(new_CutCoeffs);


    }


}