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
data = pd.read_csv('../../../../data/redx2009_full.csv', delimiter='\t', header=0, comment='#')
print(data)
sNN = np.array(data.cme)
qsq = np.array(data.q2)
x   = np.array(data.x)
dat = np.array(data.redx)
err = np.array(data.err)

# theory
def mod(x, x0, la, ga, ec, si):
    dis.set_var(x0, la, ga, ec, si, 'mve')
    return [dis.reduced_x(i, np.power(x0/i, la), 319.) for i in x]

# xx = [x0, lamb, gamma, ec, sigma]
def chi_min(xx):

    dis.set_var(xx[0], xx[1], xx[2], xx[3], xx[4], 'mve')
    the = np.array([dis.reduced_x(x[i], qsq[i], sNN[i]) for i in range(len(dat))])

    return np.sum(np.power(the - dat, 2))

# run parameters
run = 1
alg = 'pow' # fitting algorithm: 'ls', 'pow', 'shgo'

# initial guess
x0  = 1e-4
la  = 0.2
ga  = 1.
ec  = 10.
si  = 16.

# bounds order: x0, lambda, sigma/2
bounds = [(1e-5, 1e-3), (0., 0.8), (0.5, 1.5), (0., 20.), (0., 20.)]

# run minimzation
t1 = time.time()

if alg=='ls':
    print('running scipy.curve_fit...')
    popt, pcov = curve_fit(mod, x, dat, p0=[x0, la, ga, ec, si])
    print(popt)
    print(pcov)
elif alg=='pow':
    print('running scipy.optimize.minimize with Powell method...')
    res = minimize(chi_min, [x0, la, ga, ec, si], method='Powell', bounds=bounds)
    print(res)
elif alg=='shgo':
    print('running scipy.optimize.shgo...')
    # default sampling method: 'simplical', theoretical guarantee of global minimum
    # alternate methods: 'halton', 'sobol' are faster but loss of guaranteed onvergence
    res = shgo(chi_min, bounds, sampling_method='sobol')
    print(res)

t2 = time.time()

with open('out' + str(run) + '_' + alg + '.csv', 'w') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(['# run ' + str(run), 'fitting algorithm: ' + alg])
    writer.writerow(['x0', 'lamb', 'gamma', 'sigma'])
    writer.writerow([x0, la, ga, si])
    
    # write results
    if alg=='ls':
        writer.writerow(popt)
        writer.writerow(pcov)
    elif alg=='pow':
        writer.writerow(res.x)
    elif alg=='shgo':
        writer.writerow(res)
        # writer.writerow(result.x)
        # writer.writerow(['# chi2: ', results.fun])
        # writer.writerow(['# identified local minima '])
        # writer.writerow(result.xl)
        # writer.writerow(result.funl)
        # writer.writerow(['# number of iterations ', result.nit])
    # write error/confidence

    # write error/confidence
    writer.writerow(['# time: ' + str((t2-t1)/3600) + ' hours'])

