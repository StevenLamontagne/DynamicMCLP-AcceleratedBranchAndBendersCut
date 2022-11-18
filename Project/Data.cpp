#include "Data.h"



using namespace std;

void Data::load(string i, bool verbose)
{   
    if (verbose) {
        cout << "File loading starting \n";
    }
	std::ifstream f(i, ifstream::in);

    try
    {
        f >> params;
        T = (int) params["T"];
        M = (int) params["M"];
        N = (int) params["N"];
        for (long unsigned int r = 0; r < params["R"].size(); r++) {
            R.push_back((int) params["R"][r]);
        }
        for (long unsigned int k = 0; k < params["Mj"].size(); k++) {
            Mj.push_back((int) params["Mj"][k]);
        }


        if (verbose) {
            cout << "File loading success \n";
            cout << "Covering creation started \n";
        }
        try {
            create_covering(verbose);
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

double Data::SolutionQuality(vector<vector<int>> Sol)
{
    double total = 0;
    for (int t = 0; t < T; t++) {
        for (int i = 0; i < N; i++) {
            double weight = (double)params["Ni"][t][i] / (double)R[i];
            double subtotal = 0;
            for (int r = 0; r < R[i]; r++) {
                bool covered = home[t][i][r];
                for (int j = 0; j < M; j++) {
                    for (int k = 0; k < Mj[j]; k++) {
                        if ((Sol[t][j] >= k) && a[t][i][r][j][k]) {
                            covered = true;
                            break;
                        }
                    }
                    if (covered) { break; }
                }

                if (covered) { subtotal += weight; }

            }
            total += subtotal;
        }
    }
    return total;
}


void Data::create_covering(bool verbose)
{
    a.clear();
    home.clear();
    P.clear();
    cover_triplet.clear();

    //Calculate home coverage and station coverage
    for (int t = 0; t < T;t++) {
        vector<vector<vector<vector<bool>>>> a1;
        vector<vector<bool>> home1;
        for (int i = 0; i < N; i++) {
            vector<vector<vector<bool>>> a2;
            vector<bool> home2;
            for (int r = 0; r < R[i]; r++) {
                vector<vector<bool>> a3;
                double optout = params["d0"][t][0][i][r];
                if (params["d0"][t][1][i][r] == nullptr) {
                    home2.push_back(0);
                }
                else {
                    double u_home = params["d0"][t][1][i][r];
                    if (u_home >= optout) {
                        home2.push_back(1);
                    }
                    else {
                        home2.push_back(0);
                    }
                }
                for (int j = 0; j < M; j++) {
                    vector<bool> a4;
                    a4.push_back(0);
                    if (params["d1"][t][j][i][r] == nullptr) {
                        for (int k = 1; k < Mj[j]; k++) {
                            a4.push_back(0);
                        }
                    }
                    else {
                        double d1 = params["d1"][t][j][i][r];
                        for (int k = 1; k < Mj[j]; k++) {
                            if (params["beta"][t][j][i][k] == nullptr) {
                                a4.push_back(0);
                            }
                            else {
                                double beta = params["beta"][t][j][i][k];
                                double u = beta + d1;
                                if (u >= optout) {
                                    a4.push_back(1);
                                }
                                else {
                                    a4.push_back(0);
                                }
                            }

                        }            
                    }
                    a3.push_back(a4);
                }
                a2.push_back(a3);    
            }
            a1.push_back(a2);
            home1.push_back(home2); 
        }
        a.push_back(a1);
        home.push_back(home1);
    }
    

    //Precompute coverage and type of each triplet
    for (int t = 0; t < T; t++) {
        vector<vector<triplet>> P1;
        vector<vector<vector<pair<int, int>>>> cover1;
        for (int i = 0; i < N; i++) {
            vector<triplet> P2;
            vector<vector<pair<int, int>>> cover2;
            for (int r = 0; r < R[i]; r++) {
                vector<pair<int, int>> cover3;
                for (int j = 0; j < M; j++) {
                    for (int k = 1; k < Mj[j]; k++) {
                        if (a[t][i][r][j][k] == 1) {
                            cover3.push_back(make_pair(j, k));
                            break;
                        }
                    }
                }
         
                if (home[t][i][r] == 1) {
                    P2.push_back(triplet::Precovered);
                }
                else {
                    int s = cover3.size();
                    switch (s)
                    {
                    case 0:
                        P2.push_back(triplet::Uncoverable);
                        //cover3.push_back(make_pair(-1, -1));
                        break;
                    case 1:
                        P2.push_back(triplet::Single);
                        
                        break;
                    default:
                        P2.push_back(triplet::Multi);
                        break;
                    }
                }
                cover2.push_back(cover3);
            }
            P1.push_back(P2);
            cover1.push_back(cover2);
        }
        P.push_back(P1);
        cover_triplet.push_back(cover1);
    }

    //Calculate coverage and overlap
    for (int t = 0; t < T; t++) {
        vector<vector<vector<pair<int, int>>>> cover_station1;
        for (int j1 = 0; j1 < M; j1++) {
            vector<vector<pair<int, int>>> cover_station2;
            for (int k1 = 0; k1 < Mj[j1]; k1++) {
                vector<pair<int, int>> cover_station3;
                cover_station2.push_back(cover_station3);
            }
            cover_station1.push_back(cover_station2);
        }
        cover_station.push_back(cover_station1);

        for (int i = 0; i < N; i++) {
            for (int r = 0; r < R[i]; r++) {
                switch (P[t][i][r])
                {
                case triplet::Uncoverable:
                case triplet::Precovered: {
                    continue;
                    break;
                }
                default:
                {
                    for (pair<int, int> cover_triplet : cover_triplet[t][i][r]) {
                        int j = cover_triplet.first;
                        int k = cover_triplet.second;
                        //if (j >= 8|| k >= 4) {
                        //	cout << "High element found";
                        //}
                        cover_station[t][j][k].push_back(make_pair(i, r));
                    }
                    break;
                }
                }

            }
        }


        map<pair<int, int>, map<pair<int, int>, double>> overlap1;
        for (int j1 = 0; j1 < M; j1++) {
            int total = 0;
            vector<pair<int, int>> covered;
            for (int k1 = 1; k1 < Mj[j1]; k1++) {
                map<pair<int, int>, double> overlap2;
                total += cover_station[t][j1][k1].size();
                covered.insert(covered.end(), cover_station[t][j1][k1].begin(), cover_station[t][j1][k1].end());
                for (int j2 = 0; j2 < M; j2++) {
                    if (j1 == j2) { continue; }
                    for (int k2 = 1; k2 < Mj[j2]; k2++) {
                        int over = 0;
                        for (pair<int, int> trip : covered) {
                            int i = trip.first;
                            int r = trip.second;
                            if (a[t][i][r][j2][k2]) { over += 1; } //not very efficient, since recalculating a ton of stuff. But works for now
                        }
                        overlap2[make_pair(j2, k2)] = (double)over / (double)total;
                    }
                }
                overlap1[make_pair(j1, k1)] = overlap2;
            }
        }
        overlap.push_back(overlap1);


    }


}
