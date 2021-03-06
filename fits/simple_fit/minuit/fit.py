import numpy as np
import pandas as pd
from iminuit import Minuit
from iminuit.cost import LeastSquares
import time
import csv

import sys
sys.path.append('../')
import dis_solver as dis

print('packages loaded...')
import warnings
warnings.filterwarnings("ignore")

# redx data import
dnam = 'redx2009_full.csv'
data = pd.read_csv('../../../data/' + dnam, delimiter='\t', header=0, comment='#')
print(data)
sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory
# ver: parametrization (mv, mvg, mve, rcbk)
ver = 'mvg'
run = 5
x0_ = 6.5170809e-5
la_ = 0.265686
si_ = 9.60794
ga_ = 1.06990232
ec_ = 1.

print(ver)
print('run: ' + str(run))
print('init: ' + 
      'x0      = '  + str(x0_) + 
      ', lamb  = '  + str(la_) +
      ', sigma = '  + str(si_) +
      ', gamma = '  + str(ga_) +
      ', ec    = '  + str(ec_))
print('data file: ' + dnam)
def chi_squared(x0, lamb, sigma, gamma):

    if ver == 'mv':
        gamma = 1.
        ec    = 1.
    elif ver == 'mvg':
        ec    = 1.

    dis.set_var(x0, lamb, gamma, ec, sigma, ver)
    print('x0=' + str(x0) + ', lamb=' + str(lamb) + ', gamma=' + str(gamma)+ ' ec=' + str(ec)+', sigma=' + str(sigma) + ', bkver: ' + ver)
    print('vars set')
    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i])
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)
    return res
t1 = time.time()
# chi_squared.errordef = Minuit.LEAST_SQUARES
print('making instance of Minuit class')

if ver == 'mv':
    m = Minuit(chi_squared, x0=x0_, lamb=la_, sigma=si_)
elif ver =='mvg':
    m = Minuit(chi_squared, x0=x0_, lamb=la_, sigma=si_, gamma=ga_)
elif ver == 'mve':
    m = Minuit(chi_squared, x0=x0_, lamb=la_, sigma=si_, gamma=ga_, ec=ec_)

print('calling simplex method')
m.simplex()
print('simplex method complete')
t2 = time.time()
print('total fit run time: ' + str((t2 - t1)/3600) + ' hours')

print(m.values)  # prints fitted values
print(m.errors)  # prints errors
print(repr(m.fmin))

with open(ver + '_out2.txt', 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(m.values)
    writer.writerow(m.errors)
    writer.writerow(repr(m.fmin))
