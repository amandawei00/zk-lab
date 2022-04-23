# from FBT import FBT
import numpy as np
import pandas as pd
from scipy.interpolate import interp2d
from scipy.integrate import quad
from scipy.special import j0
import time
import cmath


class N:

    def __init__(self, bk, vers=''):
        self.n_ = 400
        self.x0 = 0.01
        self.xr1 = np.log(3.e-6)
        self.xr2 = np.log(60.e0)  # limit of integration in fourier transform calculation
        tol = 1.e-8
        self.width = 100
        self.pointsy = 600
        self.pointsr = 600


        # read BK solution (accepts file or pandas dataframe) to pandas dataframe
        if isinstance(bk, str):
            self.df = pd.read_csv(bk, sep="\t", comment='#', header=2)
        else:
            self.df = bk
        self.df.columns = ['y', 'vr', 'vfr']

        # converting dataframe element types
        self.df["y"] = (self.df["y"].astype('float32')).round(decimals=1)
        self.df["vr"] = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')


        self.r = np.concatenate(np.array(self.df.loc[self.df['y'] == 0.][['vr']]))  # r values for N interp
        self.y = np.unique(np.array(self.df[['y']]))

        self.z = [self.df.loc[(self.df['y'] == yy) & ((self.df['vr'] < (rr + tol)) & (self.df['vr'] > (rr - tol)))][['vfr']].iloc[0]['vfr'] for rr in self.r for yy in self.y]
        t1 = time.time()
        self.n = interp2d(self.r, self.y, self.z, kind='cubic')

        if vers != 'dis':
            t2 = time.time()
            self.k        = np.linspace(0., 8., 200)
            self.x        = [self.x0/np.exp(j) for j in self.y]
            t3 = time.time()
            self.nff_grid = [[self.udgf(kk, xx) for kk in self.k] for xx in self.x]
            t4 = time.time()
            self.nfa_grid = [[self.udga(kk, xx) for kk in self.k] for xx in self.x]
            t5 = time.time()
  
            print(self.r.shape)
            print(self.y.shape)
            print(np.array(self.z).shape)

            print(self.k.shape)
            print(np.array(self.x).shape)
            print(np.array(self.nff_grid).shape)
            # interpolate
            self.nff = interp2d(self.k, self.x, self.nff_grid, kind='cubic')
            t7 = time.time()
            self.nfa = interp2d(self.k, self.x, self.nfa_grid, kind='cubic')
            t8 = time.time()

            print('Grid of FT of N(F) took ' + str((t4-t3)/60) + ' minutes to generate')
            print('Grid of FT of N(A) took ' + str((t5-t4)/60) + ' mintues to generate')
            print('N(r,Y) took ' + str((t2-t1)/60) + ' minutes to interpolate')
            print('FT of N(F) took ' + str((t7-t5)/60) + ' minutes to interpolate')
            print('FT of N(A) took ' + str((t8-t7)/60) + ' minutes to interpolate')
       


    def get_n():
        return self.n

    # finds location of val in grid. if val is between two values of grid, lower value is returned
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

    def udgf(self, k, x):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n(r_, y_)) * j0(k * r_) * r_ 
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    def udga(self, k, x):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n_adj(r_, y_)) * j0(k * r_) * r_
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    # def bessel(self, x, alpha):
    #     f = lambda t: np.cos(alpha * t - x * np.sin(t))
    #     return (1 / np.pi) * quad(f, 0, np.pi, epsabs=0.0, epsrel=0.05)[0]

# end of class
