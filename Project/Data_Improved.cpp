#include "Data_Improved.h"



void Data_Improved::load(string path, string file, bool verbose)
{
    if (verbose) {
        cout << "File loading starting \n";
    }
    std::ifstream f(path + "Shared.json", ifstream::in);

    try
    {
        f >> params;
        T = (int)params["T"];
        

        if (verbose) {
            cout << "File loading success \n";
            cout << "Covering creation started \n";
        }
        try {
            create_covering(path, file, verbose);
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


    f.close();


}

double Data_Improved::SolutionQuality(vector<vector<int>> Sol)
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




void Data_Improved::create_covering(string path, string file, bool verbose)
{
    a.clear();
    Ps.clear();
    Precovered.clear();

    if (verbose) {
        cout << "Loading coverage file \n";
    }


    json temp_coverage;
    std::ifstream f2(path + file, ifstream::in);
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

        //Eigen::VectorXd temp_weight;
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

    //    map<pair<int, int>, map<pair<int, int>, double>> overlap1;
    //    for (int j1 = 0; j1 < M; j1++) {
    //        int total = 0;
    //        vector<pair<int, int>> covered;
    //        for (int k1 = 1; k1 < Mj[j1]; k1++) {
    //            map<pair<int, int>, double> overlap2;
    //            total += cover_station[t][j1][k1].size();
    //            covered.insert(covered.end(), cover_station[t][j1][k1].begin(), cover_station[t][j1][k1].end());
    //            for (int j2 = 0; j2 < M; j2++) {
    //                if (j1 == j2) { continue; }
    //                for (int k2 = 1; k2 < Mj[j2]; k2++) {
    //                    int over = 0;
    //                    for (pair<int, int> trip : covered) {
    //                        int i = trip.first;
    //                        int r = trip.second;
    //                        if (a[t][i][r][j2][k2]) { over += 1; } //not very efficient, since recalculating a ton of stuff. But works for now
    //                    }
    //                    overlap2[make_pair(j2, k2)] = (double)over / (double)total;
    //                }
    //            }
    //            overlap1[make_pair(j1, k1)] = overlap2;
    //        }
    //    }
    //    overlap.push_back(overlap1);
    //}


}