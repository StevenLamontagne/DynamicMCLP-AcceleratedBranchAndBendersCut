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
    cover.clear();

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
        cout << "Beginning t = " << t << endl;
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
        cout << "Completing t = " << t << endl;
        P.push_back(P1);
        cover.push_back(cover1);
        cout << "cover1 pushed" << endl;
        cout << "\n" << endl;
    }
    cout << "Triplet types computed" << endl;

}
