import pandas as pd
import numpy as np
prefix = '../redx-2009/Table'
suffix = '.csv'
df     = []
cols   = ['x', 'y', 'cme', 'sig', 'e1', 'e2', 'e3' ,'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'tot+', 'tot-']
q2     = [0.045, 0.065, 0.085, 0.1, 0.15, 0.2, 0.25, 0.35, 0.4, 0.5, 0.65, 0.85, 1.2, 1.5, 2., 2.7, 3.5, 4.5, 6.5, 8.5, 10., 12., 15., 18., 22., 27., 35., 45.]
def strip(df, col, match):
    return df.drop(df[df.x == match].index)

for i in range(1, 29):
    fname = prefix + str(i) + suffix
    panda = pd.read_csv(fname, comment='#', header=0, delimiter=',', index_col=None)

    # beginning of f2 table:
    loc = np.where(panda['X'].to_numpy() == 'X')[0][0]
    panda = panda[:loc]
    panda = panda.replace({'%':''}, regex=True)
    panda = panda.apply(pd.to_numeric)
    for j in range(len(panda.to_numpy())):
        df.append(panda.to_numpy()[j])

df1 = pd.DataFrame(df)
df1.columns = cols
print(df1)
