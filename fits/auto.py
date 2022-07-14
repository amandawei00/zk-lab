import numpy as np
import pandas as pd
import csv
import sys
sys.path.append('../bk/')
sys.path.append('../dis/')

import bk_solver as bk
from bk_interpolate import N
import dis_solver as dis

# import experimental data to pandas dataframe
def data_import(fname):
    return pd.read_csv(fname, delimiter='\t', header=0, comment='#')

# run bk with fitted parameters and write to file
def run_bk(qsq0, c2, gamma, ec, fname):
    bk_ = bk.master(qsq0, c2, gamma, ec)
    bk_.to_csv(fname, sep='\t', index=False)
    return bk_

def import_bk(fname):
    return pd.read_csv(fname, delimiter='\t', header=None, comment='#', index_col=None)

# call interpolator on bk
def bk_interp(bk_, t):
    bki = N(bk_, t)
    dis.set_n(bki)

# run dis for all data in fitdata_dis.csv and write to file
def run_dis(exp, fitted_sig):
    th_ = []
    # with open(fname, 'w') as foo:
        # writer = csv.writer(foo, delimiter='\t')
        # writer.writerow(['q2', 'cme', 'x', 'redx'])

    for i in exp.index:
        q2  = exp['q2'][i]
        cme = exp['cme'][i]
        x   = exp['x'][i]

        d   = dis.reduced_x(x, q2, cme, fitted_sig)[2]
        th_.append([q2, x, cme, d])
        # writer.writerow([q2, cme, x, d])
    df = pd.DataFrame(th_, index=None)
    df.columns = ['q2', 'x', 'cme', 'redx']
    print(len(df))
    return df

# calculate chi2/dof value with three numpy array type inputs
def chi2(th, exp, exp_err):
    err = 0.01 * exp_err * exp
    chi2 = np.sum(np.power((th - exp)/err, 2))/len(th)
    print('reduced chi2 value of fit (chi2/d.o.f): ' + str(chi2))
    return chi2

# plotting routine that takes array of Q2 values

def plot(q2):
    return 0

'''
# read results from fit file
qsq0  = 0.09179473
c2    = 10.251189017
gamma = 1.041015632
sig   = 35.04884137
ec    = 20.445469648
to_file_bk = 'fit-mve_bk.csv'
to_file_dis = 'fit-mve_dis.csv'
data = data_import('../data/fitdata_dis.csv')
# bk   = run_bk(qsq0, c2, gamma, ec, to_file_bk)
# bk_interp(bk, 'dis')
# dis  = run_dis(data, to_file_dis)
dis = data_import('results/MV/fit1_dis.csv')
chi2(dis['redx'].to_numpy(), data['redx'].to_numpy(), data['err'].to_numpy())
 
'''

# data = data_import('../data/fitdata_dis.csv')
# th   = run_dis(data, 'MV1_test_chi2.csv')
print('MVg')
print('r integral bounds: (1e-6, 1e3)')
print('r prec : 1e-8')
print('z prec : 1e-6')
print('data   : redx2009_full.csv')
print('b      : amanda')

data = data_import('../data/redx2009_full.csv')
print(data)
bk_  = import_bk('../../rcbk/mvg_test.csv')
# bk_ = import_bk('../bk/results/RK4/bk_MVg.csv')
bk_interp(bk_, 'dis')
df_  = run_dis(data, 2 * 16.45)
# df_ = pd.read_csv('heikki_hera_redx.csv', header=0, delimiter='\t')
chi2(df_['redx'], data['redx'], data['err'])

