import numpy as np
import pandas as pd
import csv
import sys

import dis_solver as dis

print('calculating chi2 for 5 parameter MVe fit')
# import experimental data to pandas dataframe
def data_import(fname):
    return pd.read_csv(fname, delimiter='\t', header=0, comment='#')

# run bk with fitted parameters and write to file
'''def run_bk(qsq0, c2, gamma, ec, fname):
    bk_ = bk.master(qsq0, c2, gamma, ec)
    bk_.to_csv(fname, sep='\t', index=False)
    return bk_
'''
def import_bk(fname):
    return pd.read_csv(fname, delimiter='\t', header=0, comment='#', index_col=None)

# call interpolator on bk
'''def bk_interp(bk_, t):
    bki = N(bk_, t)
    dis.set_n(bki)
'''
# run dis for all data in fitdata_dis.csv and write to file
def run_dis(exp, fname=''):
    th_ = []
    if fname != '':
        with open(fname, 'w') as foo:
            writer = csv.writer(foo, delimiter='\t')
            writer.writerow(['q2', 'cme', 'x', 'redx'])

            for i in exp.index:
                q2  = exp['q2'][i]
                cme = exp['cme'][i]
                x   = exp['x'][i]
            
                d   = dis.reduced_x(x, q2, cme)
                th_.append([q2, cme, x, d])
                writer.writerow([q2, cme, x, d])
    else:
        for i in range(len(exp)):
            q2  = exp['q2'][i]
            cme = exp['cme'][i]
            x   = exp['x'][i]

            d   = dis.reduced_x(x, q2, cme)
            th_.append([q2, cme, x, d])

    df = pd.DataFrame(th_, index=None)
    df.columns = ['q2', 'cme', 'x', 'redx']
    return df

# calculate chi2/dof value with three numpy array type inputs
def chi2(th, exp, exp_err):
    exp_err = 0.01 * exp_err * exp
    chi2 = np.sum(np.power((th - exp)/exp_err,2))/len(th)
    print(len(th))
    print('reduced chi2 value of fit (chi2/d.o.f): ' + str(chi2))
    return chi2

# plotting routine that takes array of Q2 values

def plot(q2):
    return 0


# read results from fit file
bk_ver= 'mv'
run   = 3
x0    = 1.08349998e-5
lamb  = 0.11377106
gamma = 1.
sig   = 7.2235975
ec    = 1.
print(bk_ver + str(run) + '_out.txt')
dis.set_var(x0, lamb, gamma, ec, sig, bk_ver)
data = data_import('../../data/redx2009_full.csv')
calculated_dis   = run_dis(data)
print(calculated_dis)
chi2(calculated_dis['redx'].to_numpy(), data['redx'].to_numpy(), data['err'].to_numpy())
 
'''
# data = data_import('../data/fitdata_dis.csv')
# th   = run_dis(data, 'MV1_test_chi2.csv')
data = data_import('../data/redx-2009/redx-2009-parsed.csv')
th_  = import_bk('../bk/results/bk_MVg2.csv')
bk_interp(th_, 'dis')
run_dis(data, '../dis/redx-2009_results/MVg2.csv', 16.45*2)
th   = data_import('../dis/redx-2009_results/MVg2.csv')
chi2(th['redx'], data['redx'], data['err'])
'''
