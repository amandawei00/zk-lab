import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
x0    = 0.01
lamb  = 0.2
gamma = 1.
ec    = 1.

x_ = np.logspace(-5, -2, 400)
r_ = np.logspace(-6, 2, 400)

def mv(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * r * r * q20
    t2 = 1/(0.241 * r) + np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))

def mvg(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * np.power(r * r * q20, gamma)
    t2 = 1/(0.241 * r) + np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))

def mve(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * r * r * q20
    t2 = 1/(0.241 * r) + ec * np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))


xxx = x_[2]
# grids
mv_grid  = pd.DataFrame([[mv(rr, xx) for xx in x_] for rr in r_])
mvg_grid = pd.DataFrame([[mvg(rr, xx) for xx in x_] for rr in r_])
mve_grid = pd.DataFrame([[mve(rr, xx) for xx in x_] for rr in r_])


mv_grid.to_csv('mv_grid.csv', sep='\t', header=False, index=False)
mvg_grid.to_csv('mvg_grid.csv', sep='\t', header=False, index=False)
mve_grid.to_csv('mve_grid.csv', sep='\t', header=False, index=False)

mv_grid_y = mv_grid[2]
mvg_grid_y = mvg_grid[2]
mve_grid_y = mve_grid[2]

# interpolated objects
a = [mv(rr, xx) for rr in r_ for xx in x_]
b = [mvg(rr, xx) for rr in r_ for xx in x_]
c = [mve(rr, xx) for rr in r_ for xx in x_]

mv_f  = interp2d(r_, x_, a, kind='cubic')
mvg_f = interp2d(r_, x_, b, kind='cubic')
mve_f = interp2d(r_, x_, c, kind='cubic')


r_grid    = np.logspace(-6, 2, 500)
mv_inter  = [mv(rr, xxx) for rr in r_grid]
mvg_inter = [mvg(rr, xxx) for rr in r_grid]
mve_inter = [mve(rr, xxx) for rr in r_grid]

plt.plot(r_, mv_grid_y, label='mv grid')
plt.plot(r_, mvg_grid_y, label='mvg grid')
plt.plot(r_, mve_grid_y, label='mve grid')

plt.plot(r_grid, mv_inter, label='mv interpolated')
plt.plot(r_grid, mvg_inter, label='mvg interpolated')
plt.plot(r_grid, mve_inter, label='mve interpolated')

plt.xscale('log')
plt.legend()
plt.show()
