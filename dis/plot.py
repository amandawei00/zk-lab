import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

qsq2   = 200
sqrt_z = 319

exp_name = '../data/reduced_x-2010.csv'
th_name  = 'x_MVg1.csv'

# import experimental data--------------------------------------
df_exp = pd.read_csv(exp_name, sep='\t', header=0, comment='#')
df_exp = df_exp.loc[df_exp['qsq2'] == qsq2]

x1  = df_exp['x']
y1  = df_exp['sig']
err = df_exp['err(%)'].multiply(0.01)

# import theoretical solutions----------------------------------
df_th  = pd.read_csv(th_name, sep='\t', header=0, comment='#')
df_th  = df_th.loc[df_th['qsq2'] == qsq2]

x2 = df_th['x']
y2 = df_th['sig']
y2 = y2.multiply(27)

plt.xlim(1.e-5, 0.15)
plt.ylim(0., 2.0)

plt.xscale('log')
plt.yscale('linear')

plt.xlabel('x', fontsize=16)
plt.ylabel('reduced cross section', fontsize=16)
plt.title('Q^2 = ' + str(qsq2) + ' GeV^2', fontsize=16)

plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
plt.errorbar(x1, y1, yerr=err, linestyle='', color='blue', capsize=2)
plt.plot(x2, y2, marker='v', linestyle='--', color='magenta', label='calculated')

plt.show()



