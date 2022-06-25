import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# script to plot experimental and theoretical dis results

order  = 'RK4'
obs    = 'redx'
qsq2   = 27.
sqrt_z = 319.
bk_ver = 'mve'

exp_name = '../../data/redx-2009-parsed.csv'
th_name  = 'mve_dis.csv'

# import experimental data--------------------------------------
df_exp = pd.read_csv(exp_name, sep='\t', header=0, comment='#')
df_exp = df_exp.loc[(df_exp['q2'] == qsq2) & (df_exp['cme'] == sqrt_z)]

x1  = df_exp['x']
y1  = df_exp[obs]
# err = df_exp['err(%)'].multiply(0.01)
err = df_exp['err'].multiply(0.01)

# import theoretical solutions----------------------------------
df_th  = pd.read_csv(th_name, sep='\t', header=0, comment='#')
df_th  = df_th.loc[(df_th['q2'] == qsq2)]

x2 = df_th['x']
f2 = df_th[obs]

plt.xlim(1.e-4, 1.e-2)
plt.ylim(0., 2.)

plt.xscale('log')
plt.yscale('linear')

plt.xlabel('x', fontsize=16)
plt.ylabel(obs, fontsize=16)
plt.title('Q^2 = ' + str(qsq2) + ' GeV^2', fontsize=16)

plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
plt.errorbar(x1, y1, yerr=err, linestyle='', color='green', capsize=2)
plt.plot(x2, f2, marker='v', linestyle='--', color='magenta', label=bk_ver, markersize=2.)
plt.legend()
plt.show()



