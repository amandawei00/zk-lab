import numpy as np
import pandas as pd
import csv
from iminuit import Minuit
from iminuit.cost import LeastSquares

import sys
sys.apth.append('../dis')
sys.path.append('../bk')
sys.path.append('../data')
import bk_solver as bk
import dis_solver as dis

import warnings
warnings.filterwarnings("ignore")

# theory

alpha = 1/137

# data import
def import_f2():

data = pd.read_csv("f2_dat.csv.csv")
data.columns = ['qsq2', 'x', 'f2', 'f2_staterr', 'f2_syserr', 'f2_experr']

qsq2 = np.array(data.qsq2)
x = np.array(data.x)
f2_dat = np.array(data.f2)
f2_err = np.array(data.f2_experr)

'''parameters:
   1. qsq20: initial saturation scale
   2. c    :
   3. gamma: anomalous dimension thing
   4. sigma: normalization factor
'''
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
m.simplex()
