import numpy as np
import pandas as pd
from scipy.optimize import shgo
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import time
import csv

import sys
sys.path.append('../../')
import dis_solver as dis

print('packages loaded...')
import warnings
warnings.filterwarnings("ignore")

# redx data import
data = pd.read_csv('../../../../data/redx-2009-parsed.csv', delimiter='\t', header=0, comment='#')
print(data)
sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory
def chi_squared(x, x0, la, ga, ec, si):
    dis.set_var(x0, la, ga, ec, si, 'mve')
    return [dis.reduced_x(i, np.power(x0/i, la), 319.) for i in x]

# xx = [x0, lamb, gamma, ec, sigma]
def chi_min(xx):

    dis.set_var(xx[0], xx[1], xx[2], xx[3], xx[4], 'mve')

    print('vars set')
    res = 0
    for i in range(len(data)):
        theory = dis.reduced_x(x[i], qsq[i], sNN[i])
        exp    = dat[i]
        er1    = err[i]
        res    += (theory - exp) * (theory - exp) / (er1 * er1)
    return res

# run parameters
run = 1
alg = 'ls' # fitting algorithm: 'ls', pow', 'shgo'

# initial guess
x0  = 1e-4
la  = 0.5
ga  = 1.
ec  = 1.
si  = 10.

# run minimization
t1 = time.time()

if alg=='ls':
    popt, pcov = curve_fit(chi_squared, x, dat, p0=[x0, la, ga, ec, si])
elif alg=='pow':
elif alg=='shgo':
    

t2 = time.time()

with open('out' + str(run) + '_' + alg + '.csv', 'w') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(['# run ' + str(run), 'fitting algorithm: ' + alg])
    writer.writerow(['x0', 'lamb', 'gamma', 'ec', 'sigma'])
    writer.writeorw([x0, la, ga, ec, si])
    # write results
    # write error/confidence
    writer.writerow(['# time: ' + str(t2-t1/3600) + ' hours')
