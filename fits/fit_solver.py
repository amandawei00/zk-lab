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
data = pd.read_csv('../data/redx2009_full.csv', delimiter='\t', header=0, index_col=None, comment='#')

print(data)
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

def chi_squared(qsq0, c, sigma):
    # run BK for given parameters qsq2, c, sigma, and gamma
    # load dataframe directly without writing to file?
    # write to file so future runs can be avoided?

    # print('set parameters: qsq0 = ' + str(qsq0) + ', c = ' + str(c) + ', g = ' + str(gamma) + ', sig = ' + str(sigma))
    print('set parameters: qsq0 = '  + str(qsq0) + ', c = ' + str(c) + ', sig = ' + str(sigma))
    bk_df = bk.master(qsq0, c, 1., 1., order='RK2')
    print('bk solution done...')
    bk_f  = N(bk_df) # why interpolate now? interpolation happens in dis
        
    # set n for dis, pp-pA
    dis.set_n(bk_f)
    print('bk interpolation done... calculating residuals')

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i], 2 * sigma)[2]
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)
    return res
t1 = time.time()
# chi_squared.errordef = Minuit.LEAST_SQUARES
print('making instance of Minuit class')
m = Minuit(chi_squared, qsq0=0.1, c=10, sigma=10)
print('calling simplex method')
m.simplex()
print('simplex method complete')
t2 = time.time()
print('total fit run time: ' + str((t2 - t1)/3600) + ' hours')
print(m.values)  # prints fitted values
print(m.errors)  # prints errors
# print(m.fval/(len(dat) - len(m.values))) # prints goodness of fit (reduced_chi2)
print(repr(m.fmin))
