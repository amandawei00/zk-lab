import numpy as np
import pandas as pd
import csv
import sys
import toydis as dis


print('calculating chi2 for 5 parameter MVe fit')
# import experimental data to pandas dataframe
def data_import(fname):
    return pd.read_csv(fname, delimiter='\t', header=0, comment='#', index_col=0)

def n(r, x, x0, lamb, gamma=1., ec=1.):
    q2 = np.power(x0/x, gamma)
    exp = -0.25 * r * r * q2 * np.log(1/(r * 0.241) + np.exp(1) * ec)
    return 1 - np.exp(exp)

# run dis for all data in fitdata_dis.csv and write to file
def run_dis(exp, fname):
    th_ = [] 
    with open(fname, 'w') as foo:
        writer = csv.writer(foo, delimiter='\t')
        writer.writerow(['q2', 'cme', 'x', 'redx'])

        for i in exp.index:
            q2  = exp['q2'][i]
            cme = exp['cme'][i]
            x   = exp['x'][i]

            d   = dis.reduced_x(x, q2, cme, sig)[2]
            th_.append(d)
            writer.writerow([q2, cme, x, d])
    return pd.DataFrame(th_, index=None)

# calculate chi2/dof value with three numpy array type inputs
def chi2(th, exp, exp_err):
    exp_err = 0.01 * exp_err * exp
    chi2 = np.sum(np.power((th - exp)/exp_err, 2))/len(th)
    print('reduced chi2 value of fit (chi2/d.o.f): ' + str(chi2))
    return chi2

# plotting routine that takes array of Q2 values

def plot(q2):
    return 0

# read results from fit 
x0    = 
lamb  = 0.41548657
gamma = 0.28681923
sig   = 

data = data_import('../data/fitdata_dis.csv')
a = lambda r, x: n(r, x, x0, lamb, gamma)
dis.set_n(a)
th = run_dis(data, 'toyfit_MV.csv')
th = data_import('toyfit_MV.csv')
chi2(th['redx'], data['redx'], data['err'])

