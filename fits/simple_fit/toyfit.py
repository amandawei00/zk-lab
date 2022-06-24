import numpy as np
import pandas as pd
import csv
from iminuit import Minuit
from iminuit.cost import LeastSquares
import time

import sys
sys.path.append('../data/')
sys.path.append('../bk/')
sys.path.append('../dis/')

import bk_solver as bk
import dis_solver as dis
from bk_interpolate import N

print('packages loaded...')
import warnings
warnings.filterwarnings("ignore")

print('fitting 5 parameters (qsq0, gamma, c, ec, sig) with no fixed parameters')
# parameters
alpha = 1/137

# data import
# f2
'''data = pd.read_csv("toydata.csv", delimiter=',', header=0, comment='#')
data.columns = ['qsq2', 'x', 'f2', 'staterr', 'syserr', 'toterr']
qsq = np.array(data.qsq2)
x   = np.array(data.x)
dat = np.array(data.f2)
err = np.array(data.toterr)'''

# reduced-x
data = pd.read_csv('../data/fitdata_dis.csv', delimiter='\t', header=0, index_col=0, comment='#')

sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory    
'''
parameters:
   1. qsq20: initial saturation scale
   2. c    :
   3. gamma: anomalous dimension thing
   4. sigma: normalization factor
'''

# parm: parametrization (mv, mvg, mve, rcbk)
def chi_squared(x0, lamb, gamma, ec, sigma):

    dis.set_var(x0, lamb, gamma, ec, sigma, 'mv')

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], sNN[i])[2]
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)
    return res
t1 = time.time()
chi_squared.errordef = Minuit.LEAST_SQUARES
print('making instance of Minuit class')
m = Minuit(chi_squared, qsq0=0.1, c=10, gamma=1., ec=20., sigma=36)
print('calling simplex method')
m.simplex()
print('simplex method complete')
t2 = time.time()
print('total fit run time: ' + str((t2 - t1)/3600) + ' hours')
print(m.values)  # prints fitted values
print(m.errors)  # prints errors
# print(m.fval/(len(dat) - len(m.values))) # prints goodness of fit (reduced_chi2)
print(repr(m.fmin))