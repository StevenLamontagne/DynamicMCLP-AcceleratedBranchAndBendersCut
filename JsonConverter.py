
import numpy as np
import os 
import json
import pickle

for set in ["Simple"]:
    for test in range(20):
        fp = "/local_1/outer/lamste/Data/Precomputed/{}/MaximumCover/MC{}".format(set, test)
        data = pickle.load(open(fp +".pickle", "rb"))
        T = data["T"]
        M = data["M"]
        N = data["N"]
        Mj = data["Mj"]

        for t in range(T):
            for j in range(M):
                for i in range(N):
                    #d1
                    for r in range(data["R"][i]):
                        if np.isnan(data["d1"][t][j][i][r]):
                            data["d1"][t][j][i][r] = None
                    #beta
                    for k in range(Mj[j]):
                        if np.isnan(data["beta"][t][j][i][k]):
                            data["beta"][t][j][i][k] = None
            #d0
            for i in range(N):
                for r in range(data["R"][i]):
                    if np.isnan(data["d0"][t][1][i][r]):
                        data["d0"][t][1][i][r] = None

        json.dump(data, open(fp +".json", "w"), allow_nan=False)




