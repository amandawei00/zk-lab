import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

qsq2   = 1.5
sqrt_z = 319.

# exp_name = '../data/redx-2009-parsed.csv'
# exp_name = '../data/redx-2009/redx-2009-parsed.csv'
# exp_name = '../data/reduced_x-2010.csv'  # good fit with x_MV1.csv
# th_name  = '../fits/results/MV/fit1_dis.csv
exp_name = '../data/redx-2009/redx-2009-parsed.csv'
# th_name  = '../dis/redx-2009_results/MV2.csv'
# exp_name   = '../data/redx-2015/fitdata_dis.csv'
th_name  = '../dis/results/redx/RK4/MV.csv'
# th_name = '../dis/results/redx/RK4/MVe_equal_masses.csv'
# import experimental data--------------------------------------
df_exp = pd.read_csv(exp_name, sep='\t', header=0, comment='#')
df_exp = df_exp.loc[(df_exp['q2'] == qsq2) & (df_exp['cme'] == sqrt_z)]

x1  = df_exp['x']
y1  = df_exp['redx']
# err = df_exp['err(%)'].multiply(0.01)
err = df_exp['err'].multiply(0.01)

# import theoretical solutions----------------------------------
df_th  = pd.read_csv(th_name, sep='\t', header=0, comment='#')
df_th  = df_th.loc[(df_th['q2'] == qsq2) & (df_th['cme'] == sqrt_z)]

x2 = df_th['x']
f2 = df_th['redx']

plt.xlim(1.e-5, 1.e-2)
plt.ylim(0., 2.)

plt.xscale('log')
plt.yscale('linear')

plt.xlabel('x', fontsize=16)
plt.ylabel('reduced cross section', fontsize=16)
plt.title('Q^2 = ' + str(qsq2) + ' GeV^2', fontsize=16)

plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
plt.errorbar(x1, y1, yerr=err, linestyle='', color='blue', capsize=2)
plt.plot(x2, f2, marker='v', linestyle='--', color='magenta', label='MVg3 results with Qs^2', markersize=2.)
plt.legend()
plt.show()



