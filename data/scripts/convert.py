import numpy as np
import pandas as pd
import csv

# df = pd.read_csv('f2-combined.csv', sep=',', header=0)
# df.to_csv('f2-combined.csv', sep='\t')

err = ['e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'e15', 'e16', 'e17', 'e18', 'e19', 'e20']
df = pd.read_csv('reduced_x.csv', sep='\t', header=0, index_col=0)
df['toterr'] = df['e1'] * df['e1']

for i in range(len(err)):
    df['toterr'] += df[err[i]] * df[err[i]]

df['toterr'] = df['toterr'].pow(1./2)
print(df)
