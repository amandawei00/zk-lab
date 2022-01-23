import numpy as np
import pandas as pd
import csv
# from iminuit import Minuit
# from iminuit.cost import LeastSquares

import sys
sys.path.append('../dis')
sys.path.append('../bk')
sys.path.append('../data')
import bk_solver as bk
import dis_solver as dis

import warnings
warnings.filterwarnings("ignore")
print("import done ")
# theory

alpha = 1/137

# data import
error_name = ['e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'e15', 'e16', 'e17', 'e18', 'e19', 'e20']
data = pd.read_csv("../data/reduced_x.csv", delimiter='\t')
data.columns = ['index','qsq2', 'x', 'sig'] + error_name

# total error calculation
data['err(tot)'] = pd.DataFrame(np.zeros(len(data.index)))
for i in error_name:
    data['err(tot)'] += data[i] * data[i]
data['err(tot)'] = data['err(tot)'].apply(np.sqrt)

# print(data)
'''
qsq2 = np.array(data.qsq2)
x = np.array(data.x)

parameters:
   1. qsq20: initial saturation scale
   2. c    :
   3. gamma: anomalous dimension thing
   4. sigma: normalization factor

def chi_squared(qsq20, c, gamma, sigma):
    # run BK for given parameters qsq2, c, and gamma
    res = 0
    for i in range(len(f2_dat)):
        theory = f2(x[i], qsq2[i])
        exp = f2_dat[i]
        err = f2_err[i]
        res += (theory - exp) * (theory - exp) / (err * err)
    return res

chi_squared.errordef = Minuit.LEAST_SQUARES
m = Minuit(chi_squared, qsq20=10.08, c=1, gamma=1)
m.simplex()'''
