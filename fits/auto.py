import numpy as np
import pandas as pd
import csv
import sys
sys.path.append('../bk/')
sys.path.append('../dis/')

import bk_solver as bk
from bk_interpolate import N
import dis_solver as dis

print('calculating chi2 for 5 parameter MVe fit')
# import experimental data to pandas dataframe
exp    = pd.read_csv('../data/fitdata_dis.csv', delimiter='\t', header=0, comment='#', index_col=0)
# exp_q2 = np.unique(np.array(exp[['q2']]))

# read results from fit results file

qsq0  = 0.09179473
c2    = 10.251189017
gamma = 1.041015632
sig   = 35.04884137
ec    = 20.445469648

# run bk with fitted parameters

bk_ = bk.master(qsq0, c2, gamma, ec)
bk_.to_csv('fit-mve_bk.csv', sep='\t', index=False)
# call interpolator on bk

bki = N(bk_, 'dis')
dis.set_n(bki)

# run dis for all data in fitdata_dis.csv and write to file
th_ = []
with open('fit-mve_dis.csv', 'w') as foo:
    writer = csv.writer(foo, delimiter='\t')
    writer.writerow(['q2', 'cme', 'x', 'redx'])

    for i in exp.index:
        q2  = exp['q2'][i]
        cme = exp['cme'][i]
        x   = exp['x'][i]

        d   = dis.reduced_x(x, q2, cme, sig)[2]
        th_.append(d)
        writer.writerow([q2, cme, x, d])

# calculate chi2/dof value

exp_ = exp['redx'].to_numpy()
chi2 = np.sum(np.power(exp_ - np.array(th_),2)/np.array(th_))/len(th_)
print('reduced chi2 value of fit (chi2/d.o.f): ' + str(chi2))

# plotting routine that takes array of Q2 values



