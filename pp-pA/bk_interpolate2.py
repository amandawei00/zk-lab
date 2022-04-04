from scipy.special import j0
from scipy.fft import fft, fft2, fftfreq, fftshift
import numpy as np
import pandas as pd
from scipy.interpolate import interp2d
from scipy.integrate import quad
import cmath
import matplotlib.pyplot as plt
import time

class N:

    def __init__(self, filename):
        self.n_      = 400
        self.x0      = 0.01
        self.xr1     = np.log(3.e-6)
        self.xr2     = np.log(60.e0)  # limit of integration in fourier transform calculation
        self.width   = 100
        self.pointsy = 600
        self.pointsr = 600
        tol          = 1e-8

        # read results.csv file from BK solution to pandas dataframe
        self.df = pd.read_csv(filename, sep="\t") 
        self.df.columns = ['y', 'vr', 'vfr']

        # converting dataframe element types
        self.df["y"]   = (self.df["y"].astype('float32')).round(decimals=1)
        self.df["vr"]  = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')

        self.r = np.concatenate(np.array(self.df.loc[self.df['y'] == 0.][['vr']]))  # r values for N interp
        self.y = np.unique(np.array(self.df[['y']]))

        self.z      = [self.df.loc[(self.df['y'] == yy) & ((self.df['vr'] < (rr + tol)) & (self.df['vr'] > (rr - tol)))][['vfr']].iloc[0]['vfr'] for rr in self.r for yy in self.y]
        self.n      = interp2d(self.r,    self.y,    self.z,   kind='cubic')

        self.k      = np.linspace(2, 9, 200)
        t1 = time.time()
        self.zft_f  = [self.udgf(kk, yy) for kk in self.k for yy in self.y]
        print('udg_f grid calculated, took ' + str((time.time()-t1)/60) + ' minutes')
        t2 = time.time()
        self.zft_a  = [self.udga(kk, yy) for kk in self.k for yy in self.y] 
        print('udg_a grid calculated, took ' + str((time.time()-t2)/60) + ' minutes')

        self.uf  = interp2d(self.k, self.y, self.zft_f, kind='cubic')
        self.ua  = interp2d(self.k, self.y, self.zft_a, kind='cubic')
        
        # self.zft  = fft(self.z)
        # self.zftr = fftfreq(len(self.r))

        # self.zft  = fftshift(self.zft)
        # self.zftr = fftshift(self.zftr)
        # self.nft  = interp2d(self.zftr, self.y, self.zft, kind='cubic')

    # if val is between two elements in grid, find_index returns index of lower element
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

    # returns interpolated value of N(r,Y)
    def n_adj(self, r, y):
        return 2 * self.n(r, y) - np.power(self.n(r, y), 2)

    def udgf(self, k, y):
        integrand = lambda r_: (1 - self.n(r_, y)) * self.bessel(k * r_, 0) * r_
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    def udga(self, k, y):
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.n_adj(r_, y_)) * self.bessel(k * r_, 0) * r_
        return 2 * np.pi * quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=0.0, epsrel=0.05)[0]

    def udg_f(self, x, k):
        y_ = np.log(self.x0 / x)
        return self.uf(k, y_)

    def udg_a(self, x, k):
        y_ = np.log(self.x0 / x)
        return self.ua(k, y_)

    def bessel(self, x, alpha):
        f = lambda t: np.cos(alpha * t - x * np.sin(t))
        return (1 / np.pi) * quad(f, 0., np.pi, epsabs=0.0, epsrel=0.05)[0]
# end of class

# if __name__ == '__main__':
#     n = N('../bk/results/bk_MV1.csv')
    
