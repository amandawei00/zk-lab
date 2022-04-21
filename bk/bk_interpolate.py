# from FBT import FBT
import numpy as np
import pandas as pd
import scipy.interpolate as interp
from scipy.integrate import quad
import time
import cmath


class N:

    def __init__(self, bk):
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
            self.df = pd.read_csv(bk, sep="\t", comment='#')
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
        self.n = interp.interp2d(self.r, self.y, self.z, kind='cubic')

    # finds location of val in grid in the sense that if val is between two elements in grid, find_index will return the index of the lower element
    def get_n():
        return self.n

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
        integrand = lambda r_: (1 - self.n(r_, y_)) * self.bessel(k * r_, 0) * r_ 
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    def udga(self, x, k):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n_adj(r_, y_)) * self.bessel(k * r_, 0) * r_
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    def bessel(self, x, alpha):
        f = lambda t: np.cos(alpha * t - x * np.sin(t))
        return (1 / np.pi) * quad(f, 0, np.pi, epsabs=0.0, epsrel=0.05)[0]

# end of class
