import numpy as np
import pandas as pd
import csv
from iminuit import Minuit
from iminuit.cost import LeastSquares
import time

import sys
sys.path.append('../data/')
sys.path.append('../dis/')

import dis_solver as dis

print('packages loaded...')
import warnings
warnings.filterwarnings("ignore")

model = 'MV'
# model = 'MVg'
# model = 'MVe'

print('running fit to DIS data with ' + model + ' parametrization')

# parameters
alpha = 1/137
l_QCD = 0.241

# data import, reduced cross section DIS
data = pd.read_csv('../data/fitdata_dis.csv', delimiter='\t', header=0, index_col=0, comment='#')

sNN = np.array(data.cme)
qsq = np.array(data.qsq2)
x   = np.array(data.x)
dat = np.array(data.sig)
err = np.array(data.err)

# theory


def theory(r, x, x0, lamb, sig, gamma=1, ec=1.):
    exp = -0.24 * r 

# def chi_squared(qsq0, c, gamma, sigma):
def chi_squared(qsq0, c, gamma, ec, sigma):
    # run BK for given parameters qsq2, c, sigma, and gamma
    # load dataframe directly without writing to file?
    # write to file so future runs can be avoided?

    # print('set parameters: qsq0 = ' + str(qsq0) + ', c = ' + str(c) + ', g = ' + str(gamma) + ', sig = ' + str(sigma))
    print('set parameters: qsq0 = '  + str(qsq0) + ', c = ' + str(c) + ', sig = ' + str(sigma))
    bk_df = bk.master(qsq0, c, gamma, ec)

    print('bk solution done...')
    bk_f  = N(bk_df) # why interpolate now? interpolation happens in dis
    # set n for dis, pp-pA
    dis.set_n(bk_f)
    print('bk interpolation done... calculating residuals')

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i], sigma)[2]
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
