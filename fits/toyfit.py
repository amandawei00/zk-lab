import numpy as np
import pandas as pd
import csv
from iminuit import Minuit
from iminuit.cost import LeastSquares
import time

import sys
sys.path.append('../data/')
sys.path.append('../dis/')

import toydis as dis

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
# data = pd.read_csv('../fits/fitdata_dis.csv', delimiter='\t', header=0, index_col=0, comment='#')
data = pd.read_csv('toydata-2009.csv', delimiter='\t', header=0, index_col=None, comment='#')

sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory
def n(r, x, x0, lamb, gamma=1, ec=1.):
    q2  = np.power(x0/x, lamb)
    exp = -0.25 * np.power(r * r * q2, gamma) * np.log(1/(0.241 * r) + ec * np.exp(1))
    return 1 - np.exp(exp)

def chi2_mv(x0, lamb, sig):
    print('set parameters: x0 = ' + str(x0) + ', lamb = ' + str(lamb) +  ', sig = ' + str(sig))
    a = lambda r, x: n(r, x, x0, lamb)
    dis.set_n(a)

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i], sig)[2]
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)

    return res


def chi2_mvg(x0, lamb, gamma, sig):
    print('set parameters: x0 = ' + str(x0) + ', lamb = ' + str(lamb) +  ', gamma = ' + str(gamma) + ', sig = ' + str(sig))
    a = lambda r, x: n(r, x, x0, lamb, gamma)
    dis.set_n(a)

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i], sig)[2]
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)

    return res


def chi2_mve(x0, lamb, gamma, ec, sig):
    print('set parameters: x0 = ' + str(x0) + ', lamb = ' + str(lamb) + ', gamma = ' + str(gamma) + ', ec = ' + str(ec) + ', sig = ' + str(sig))
    a = lambda r, x: n(r, x, x0, lamb, gamma, ec)
    dis.set_n(a)

    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i], sig)[2]
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)

    return res
 
t1 = time.time()

x0_ = 1.e-4
l_  = 0.3
g_  = 1.
ec_ = 1.
s_  = 20

print('Making fit with the following initials: ')
print('x0     = ' + str(x0_))
print('lambda = ' + str(l_))
print('gamma  = ' + str(g_))
print('ec     = ' + str(ec_))

print('normalization constant (sigma) = ' + str(s_))
if model == 'MV':
    print('fitting MV parametrization (x0, lambda, sigma)')
    chi2_mv.errordef = Minuit.LEAST_SQUARES
    m = Minuit(chi2_mv, x0=x0_, lamb=l_, sig=s_)
elif model == 'MVg':
    print('fitting MVg parametrization (x0, lambda, gamma, sigma)')
    chi2_mvg.errordef = Minuit.LEAST_SQUARES
    m = Minuit(chi2_mvg, x0=x0_, lamb=l_, gamma=g_, sig=s_)
elif model == 'MVe':
    print('fitting MVe parametrization (x0, lambda, gamma, ec, sigma)')
    chi2_mve.errordef = Minuit.LEAST_SQUARES
    m = Minuit(chi2_mve, x0=x0_, lamb=l_, gamma=g_, ec=ec_, sig=s_)

print('using simplex minimization routine')
m.simplex()
t2 = time.time()
print('total fit run time: ' + str((t2 - t1)/3600) + ' hours')
print(m.values)  # prints fitted values
print(m.errors)  # prints errors
# print(m.fval/(len(dat) - len(m.values))) # prints goodness of fit (reduced_chi2)
print(repr(m.fmin))
