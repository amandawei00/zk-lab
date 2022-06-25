import numpy as np
import pandas as pd
from iminuit import Minuit
from iminuit.cost import LeastSquares
import time
import csv
import dis_solver as dis

print('packages loaded...')
import warnings
warnings.filterwarnings("ignore")

# redx data import
data = pd.read_csv('redx-2009-parsed.csv', delimiter='\t', header=0, comment='#')
print(data)
sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory
# ver: parametrization (mv, mvg, mve, rcbk)
ver = 'mve'
def chi_squared(x0, lamb, sigma, gamma, ec):

    if ver == 'mv':
        gamma = 1.
        ec    = 1.
    elif ver == 'mvg':
        ec    = 1.

    dis.set_var(x0, lamb, gamma, ec, sigma, ver)

    print('vars set')
    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i])
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)
    return res
t1 = time.time()
chi_squared.errordef = Minuit.LEAST_SQUARES
print('making instance of Minuit class')

if ver == 'mv':
    m = Minuit(chi_squared, x0=1e-4, lamb=0.5, sigma=10.)
elif ver =='mvg':
    m = Minuit(chi_squared, x0=1e-4, lamb=0.5, sigma=10., gamma=1.)
elif ver == 'mve':
    m = Minuit(chi_squared, x0=1e-4, lamb=0.5, sigma=10., gamma=1., ec=1.)

print('calling simplex method')
m.simplex()
print('simplex method complete')
t2 = time.time()
print('total fit run time: ' + str((t2 - t1)/3600) + ' hours')

print(m.values)  # prints fitted values
print(m.errors)  # prints errors
print(repr(m.fmin))

with open(ver + '_out.txt', 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(m.values)
    writer.writerow(m.errors)
    writer.writerow(repr(m.fmin))
