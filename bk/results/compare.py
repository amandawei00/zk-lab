import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dfe = pd.read_csv('bk_MV_euler.csv', sep='\t', header=2)
df  = pd.read_csv('bk_MV.csv', sep='\t', header=2)

dfe.columns = ['y', 'r', 'N(r,y)']
df.columns  = ['y', 'r', 'N(r,y)']

xe = dfe.loc[dfe['y']==0.2]['N(r,y)'].to_numpy()
xx = df.loc[df['y']==0.2]['N(r,y)'].to_numpy()

diff = (xe - xx)/xe

plt.xscale('log')
plt.plot(diff)
plt.show()
