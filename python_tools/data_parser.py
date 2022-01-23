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

def make_df(name_list):
    df_list = [parse_f2(i) for i in name_list]
    return pd.concat(df_list, ignore_index=True)

a = '../data/f2/f2_1.csv'
b = '../data/f2/f2_2.csv'
c = '../data/f2/f2_3.csv'

# make_df([a, b, c]).to_csv('f2.csv', sep='\t')
        
