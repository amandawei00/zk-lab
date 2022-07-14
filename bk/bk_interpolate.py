# from FBT import FBT
import numpy as np
import pandas as pd
from scipy.interpolate import interp2d
from scipy.integrate import quad
from scipy.special import j0
import matplotlib.pyplot as plt
import time
import cmath


class N:

    def __init__(self, bk, vers=''):
        self.x0 = 0.01
        self.r2 = 1.e2

        tol = 1.e-8
        self.pointsy = 600
        self.pointsr = 600


        # read BK solution (accepts file or pandas dataframe) to pandas dataframe
        if isinstance(bk, str):
            self.df   = pd.read_csv(bk, sep='\t', comment='#', header=None)
            self.name = bk
        else:
            self.df   = bk
            self.name = ''
        self.df.columns = ['y', 'vr', 'vfr']

        # converting dataframe element types
        # change decimals to match digits of hy
        self.df["y"] = (self.df["y"].astype('float32')).round(decimals=2)
        self.df["vr"] = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')

        # self.r = np.concatenate(np.array(self.df.loc[self.df['y'] == 0.0][['vr']]))
        self.r = np.unique(np.array(self.df[['vr']]))
        self.y = np.unique(np.array(self.df[['y']]))
        self.z = [self.df.loc[(self.df['y'] == yy) & ((self.df['vr'] < (rr + tol)) & (self.df['vr'] > (rr - tol)))][['vfr']].iloc[0]['vfr'] for rr in self.r for yy in self.y]
        t1 = time.time()

        self.nfunc = interp2d(self.r, self.y, self.z, kind='cubic')

        if vers != 'dis':
            t2 = time.time()
            self.k        = np.linspace(0., 8., 200)
            self.x        = [self.x0/np.exp(j) for j in self.y]
            t3 = time.time()
            self.nff_grid = [[self.udgf(xx, kk) for xx in self.x] for kk in self.k]
            t4 = time.time()
            self.nfa_grid = [[self.udga(xx, kk) for xx in self.x] for kk in self.k]
            t5 = time.time()
  
            print(self.r.shape)
            print(self.y.shape)
            print(np.array(self.z).shape)

            print(self.k.shape)
            print(np.array(self.x).shape)
            print(np.array(self.nff_grid).shape)
            # interpolate
            self.nff = interp2d(self.x, self.k, self.nff_grid, kind='cubic')
            t7 = time.time()
            self.nfa = interp2d(self.x, self.k, self.nfa_grid, kind='cubic')
            t8 = time.time()

            print('Grid of FT of N(F) took ' + str((t4-t3)/60) + ' minutes to generate')
            print('Grid of FT of N(A) took ' + str((t5-t4)/60) + ' mintues to generate')
            print('N(r,Y) took ' + str((t2-t1)/60) + ' minutes to interpolate')
            print('FT of N(F) took ' + str((t7-t5)/60) + ' minutes to interpolate')
            print('FT of N(A) took ' + str((t8-t7)/60) + ' minutes to interpolate')
       


    def get_n(self):
        return self.n

    def get_name(self):
        return self.name

    def n(self, r, y):
        if r > self.r2:
            return 1.
        else:
            if y < 0.0:
                print('!!!y is below lowest value in grid!!!')
            return self.nfunc(r, y)

   # finds loc of val in grid. if val is between two values of grid, returns lower value
    def find_index(self, val, grid):
        index = 0

        cont = True
        i = 0
        while cont:
            if val >= grid[i]:
                index = i
                cont = False

            i += 1

            if i >= len(grid):
                cont = False

        return index

    def n_adj(self, r, y):
        return 2 * self.n(r, y) - np.power(self.n(r, y), 2)

    def udgf(self, x, k):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n(r_, y_)) * j0(k * r_) * r_ 
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=1e-4)[0]

    def udga(self, x, k):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n_adj(r_, y_)) * j0(k * r_) * r_
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=1e-4)[0]


if __name__ == '__main__':
    y     = 5.
    x     = 0.01/np.exp(y)

    fname = 'results/RK4/bk_MV.csv'

    # grid
    df    = pd.read_csv(fname, header=None, sep='\t', comment='#')
    df.columns = ['y', 'r', 'N(r,y)']

    r     = df.loc[df['y'] == y][['r']]
    # y     = df.loc[df['y'] == y][['N(r,y)']]
    k     = np.linspace(0., 8., 200)

    # interpolated object
    bk    = N(fname)
     
    ngrid1 = [bk.udgf(x, k[i]) for i in range(len(k))]
    ngrid2 = [bk.nff(x, k[i]) for i in range(len(k))]

    plt.plot(k, ngrid1, label='grid')
    # plt.plot(k, ngrid2, label='interpolated')

    plt.xscale('log')
    plt.show()
# end of class
