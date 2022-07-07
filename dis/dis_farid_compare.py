import numpy as np
import pandas as pd
import sys
import csv

sys.path.append('../bk/')
from bk_interpolate import N
import dis_solver as dis

bk = N('../../rcbk/results/mve_test.csv', 'dis')
dis.set_n(bk)

dff = pd.read_csv('mve_dis_farid.csv', delimiter='\t', header=0)
for i in range(len(dff)):
    q2 = dff['q2'][i]
    x  = dff['x'][i]
    cme = dff['cme'][i]
    red = dff['redx'][i]

    amanda = dis.reduced_x(x, q2, cme, 16.36 * 2)[2]
    upper = red + 1e-7
    lower = red - 1e-7

    if amanda > upper or amanda < lower:
        print((red - amanda)/red)
