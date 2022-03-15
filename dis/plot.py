import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

qsq2   = 200
sqrt_z = 319

exp_name = '../data/reduced_x-2010.csv'
th_name  = 'x_MV1.csv'

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
fl = df_th['fl']
f2 = df_th['f2']
y2 = df_th['sig']

m  = 27
fl = fl.multiply(m)
f2 = f2.multiply(m)
y2 = y2.multiply(m)


plt.xlim(1.e-5, 1.e-2)
plt.ylim(0., 2.0)

plt.xscale('log')
plt.yscale('linear')

plt.xlabel('x', fontsize=16)
plt.ylabel('reduced cross section', fontsize=16)
plt.title('Q^2 = ' + str(qsq2) + ' GeV^2', fontsize=16)

plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
plt.errorbar(x1, y1, yerr=err, linestyle='', color='blue', capsize=2)
plt.plot(x2, y2, marker='v', linestyle='-', color='black', label='calculated', markersize=2.)
plt.plot(x2, fl, marker='v', linestyle='-.', color='green', label='fl', markersize=2.)
plt.plot(x2, f2, marker='v', linestyle='--', color='magenta', label='f2', markersize=2.)
plt.legend()
plt.show()



