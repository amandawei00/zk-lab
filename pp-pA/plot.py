import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

y      = 3.8
sqrt_s = 200  # [GeV]

exp_name = '../data/pp.csv'
th_name  = 'results/pp/RK2/pi0_33.csv'
# th_name  = 'pp_MV1.csv'
# import experimental data--------------------------------------
df_exp = pd.read_csv(exp_name, sep='\t', header=0, comment='#')
df_exp = df_exp.loc[df_exp['y'] == y]

x1  = df_exp['pt']
y1  = df_exp['dN']
# err = df_exp['err(%)'].multiply(0.01)

# import theoretical solutions----------------------------------
df_th  = pd.read_csv(th_name, sep='\t', header=0, comment='#')
df_th  = df_th.loc[df_th['y'] == y]

x2 = df_th['pt']
y2 = df_th['dN']

plt.xlim(0., 5.)
plt.ylim(1.e-7, 1.e2)

plt.xscale('linear')
plt.yscale('log')

plt.xlabel('x', fontsize=16)
plt.ylabel('differential cross section', fontsize=16)
plt.title('pi0, y = ' + str(y), fontsize=16)

plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
# plt.errorbar(x1, y1, yerr=err, linestyle='', color='blue', capsize=2)
plt.plot(x2, y2, marker='v', linestyle='--', color='magenta', label='calculated')

plt.show()



