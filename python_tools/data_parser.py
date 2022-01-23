import numpy as np
import pandas as pd
import csv
import sys

sys.path.append('../data/f2')
sys.path.append('../data/fl')
sys.path.append('../data/reduced-x')

# returns df stripped of anything except numbers (labels, '%',...)
def strip(df, col, match):
    return df.drop(df[df.x == match].index)

# returns dataframe of numbers so fit_solver can parse through
def parse_f2(filename):
    df = pd.read_csv(filename, delimiter=',', comment='#', header=None)
    df.columns = ['x', 'f2', 'tot(+)', 'tot(-)', 'stat(+)', 'stat(-)']
    return strip(df, 'x', 'X')

def parse_reducedx(filename):
    df = pd.read_csv(filename, delimiter=',', comment='#', header=None)
    df.columns = ['qsq2', 'x', 'sig', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'e15', 'e16', 'e17', 'e18', 'e19', 'e20'] 
    return df

def make_df(f, name_list):
    df_list = [f(i) for i in name_list]
    return pd.concat(df_list, ignore_index=True)

# filters out entries with x > x_, the small-x limit
def filter(df, x_=0.01):
    df = df.drop(df[df.x > x_].index)
    return df

f2 = '../data/f2/'
rx = '../data/reduced-x/'

a = rx + 'reduced_x1.csv'
b = rx + 'reduced_x2.csv'
c = rx + 'reduced_x3.csv'
d = rx + 'reduced_x4.csv'
e = rx + 'reduced_x5.csv'

# make_df(parse_f2, [a, b, c]).to_csv('f2.csv', sep='\t')
df = make_df(parse_reducedx, [a, b, c, d, e])
df = filter(df)
df.to_csv('reduced_x.csv', sep='\t')
