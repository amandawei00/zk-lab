import csv
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', None, 'display.max_columns', None)

df = pd.read_csv('reduced_x.csv', header=0, delimiter='\t', index_col=0)
df = df.drop(df[df.qsq2 > 50.].index)
df = df.drop(df[df.x > 0.01].index)
err_arr = [df.e1.to_numpy(), df.e2.to_numpy(), df.e3.to_numpy(), df.e4.to_numpy(), df.e5.to_numpy(), df.e6.to_numpy(), df.e7.to_numpy(), df.e8.to_numpy(), df.e9.to_numpy(), df.e10.to_numpy(), df.e11.to_numpy(), df.e12.to_numpy(), df.e13.to_numpy(), df.e14.to_numpy(), df.e15.to_numpy(), df.e16.to_numpy(), df.e17.to_numpy(), df.e18.to_numpy(), df.e19.to_numpy(), df.e20.to_numpy()]

err_arr = np.abs(np.transpose(err_arr))
max_err = [max(err_arr[i]) for i in range(len(err_arr))]
df = df.drop(['e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'e15', 'e16', 'e17', 'e18', 'e19', 'e20'], axis=1)
df.reset_index(drop=True, inplace=True)
df = df.join(pd.DataFrame(max_err))
print(df)
df.to_csv('fitdata_dis.csv', sep='\t')
