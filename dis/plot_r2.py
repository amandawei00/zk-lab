import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# script to plot experimental and theoretical dis results
obs      = 'redx'
qsq2     = 1.5
sqrt_z   = 319.

prefix   = 'results/redx/RK4/MVe_r2-'
th1      = prefix + '9e2.csv'
th2      = prefix + '95e2.csv'
th3      = prefix + '1e3.csv'
th4      = prefix + '15e3.csv'
th5      = prefix + '2e3.csv'
th6      = prefix + '3e3.csv'

th       = [th1, th2, th3, th4, th5, th6]
exp_name = '../data/redx-2009-parsed.csv'
# import experimental data--------------------------------------
df_exp = pd.read_csv(exp_name, sep='\t', header=0, comment='#')
df_exp = df_exp.loc[(df_exp['q2'] == qsq2) & (df_exp['cme'] == sqrt_z)]

x1  = df_exp['x']
y1  = df_exp[obs]
# err = df_exp['err(%)'].multiply(0.01)
err = df_exp['err'].multiply(0.01)

# import theoretical solutions----------------------------------
dfs = []
for i in range(len(th)):
    df = pd.read_csv(th[i], sep='\t', header=0, comment='#')
    dfs.append(df.loc[(df['q2'] == qsq2) & (df['cme'] == sqrt_z)])

xs = []
ys = []

for i in range(len(th)):
    xs.append(dfs[i]['x'])

for i in range(len(th)):
    ys.append(dfs[i][obs])

plt.xlim(1.e-5, 1.e-2)
plt.ylim(0., 2.)

plt.xscale('log')
plt.yscale('linear')

plt.xlabel('x', fontsize=16)
plt.ylabel(obs, fontsize=16)
plt.title('Q^2 = ' + str(qsq2) + ' GeV^2', fontsize=16)

for i in range(len(th)):
    plt.plot(xs[i], ys[i], label=th[i])

plt.legend()
plt.plot(x1, y1, marker='D', linestyle='', color='blue', label='experimental')
plt.errorbar(x1, y1, yerr=err, linestyle='', color='blue', capsize=2)

plt.legend()
plt.show()



