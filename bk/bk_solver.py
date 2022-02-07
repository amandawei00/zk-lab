import numpy as np
import csv
from multiprocessing import Pool
# import matplotlib.pyplot as plt

import time
import numba
from numba import jit
from numba import cfunc, carray
from numba.types import intc, intp, float64, voidptr
from numba.types import CPointer
from scipy import LowLevelCallable as llc
import scipy.interpolate as interpolate
from scipy.interpolate import UnivariateSpline
from scipy.integrate import nquad
from scipy.integrate import dblquad
import solver as so
import subprocess

import warnings
warnings.filterwarnings('ignore')

# variables
n = 399  # number of r points to be evaluated at each evolution step in Y
r1 = 3.e-6  # limits of r
r2 = 60.e0

xr1 = np.log(r1)
xr2 = np.log(r2)

hr = (xr2 - xr1) / n

# hy = 0.1
hy = 0.2
ymax = 10.0
y = np.arange(0.0, ymax, hy)

# Arrays for N and r in N(r), evaluated at some rapidity Y (including next step N(r,Y) in the evolution
xlr_ = [xr1 + i * hr for i in range(n + 1)]
pi_ = np.linspace(0.0, 0.5*np.pi, 20)
r_ = np.exp(xlr_)
n_ = []

# parameters
nc = 3        # number of colors
nf = 3        # number of active flavors
lamb = 0.241  # lambda QCD (default)

beta = (11 * nc - 2. * nf)/(12 * np.pi)
afr = 0.7     # frozen coupling constant (default)
rfr = (2./lamb) * np.exp(-0.5/(beta * afr))  # IR cutoff
2
c2, gamma, qs02 = 0. , 0., 0.   # fitting parameters
e = np.exp(1)
ec = 1

# initial condition
# @jit(float64(float64), nopython=True)
def mv(r):
    xlog = np.log(1/(lamb * r) + ec * e)
    xexp = np.power(qs02 * r * r, gamma) * xlog/4.0
    return 1 - np.exp(-xexp)

def evolve(xlr):

    index = xlr_.index(xlr)
    nr0 = n_[index]

    # set_vars(xlr, n(r0), xlr_, n_arr)
    so.set_vars(xlr, nr0, xlr_, n_)

    fk = llc.from_cython(so, 'f_kernel', signature='double (int, double *)')
    fs = llc.from_cython(so, 'f_split', signature='double (int, double *)')
    fc = llc.from_cython(so, 'f_combined', signature='double (int, double *)')

    Ker = dblquad(fk, xr1, xr2, 0.0, 0.5 * np.pi, epsabs=0.00, epsrel=0.05)[0]
    Spl = dblquad(fs, xr1, xr2, 0.0, 0.5 * np.pi, epsabs=0.00, epsrel=0.05)[0]
    Com = dblquad(fc, xr1, xr2, 0.0, 0.5 * np.pi, epsabs=0.00, epsrel=0.05)[0]

    k1 = Com
    k2 = k1 + (0.5 * hy * k1 * Ker) - (0.5 * hy * k1 * Spl) - (0.25 * hy * hy * k1 * k1 * Ker)
    k3 = k1 + (0.5 * hy * k2 * Ker) - (0.5 * hy * k2 * Spl) - (0.25 * hy * hy * k2 * k2 * Ker)
    k4 = k1 + (0.5 * hy * k3 * Ker) - (0.5 * hy * k3 * Spl) - (0.25 * hy * hy * k3 * k3 * Ker)

    return (1/6) * hy * (k1 + 2 * k2 + 2 * k3 + k4)

# pass fitting variables q_, c_, g_ to set variables in master.py
def master(q_, c_, g_, filename):
    global n_, qs02, c2, gamma

    # variables
    qs02 = q_
    c2 = c_
    gamma = g_

    so.set_params(c2, gamma, qs02) 

    # opening file 'results.csv' to store data from this run
    with open(filename, 'w') as csv_file:
        writer = csv.writer(csv_file, delimiter="\t")
        writer.writerow(["y", "r", "N(r,Y)"])

        # initial condition----------------------------------------------------------
        n_ = [mv(r_[i]) for i in range(len(r_))]
        #----------------------------------------------------------------------------
        # begin evolution
        for i in range(len(y)):
            y0 = y[i]
            print("y = " + str(y0))

            for j in range(len(r_)):
                print('r = ' + str(r_[j]) + ', N(r,Y) = ' + str(n_[j]))
                writer.writerow([y0, r_[j], n_[j]])
            # calculate correction and update N(r,Y) to next step in rapidity

            xk = []
            with Pool(processes=4) as pool:
                xk = pool.map(evolve, xlr_, chunksize=100)

            # xk = [evolve(xlr_[i], n_) for i in range(len(xlr_))]
            n_ = [n_[j] + xk[j] for j in range(len(n_))]

            # remove nan values from solution
            xx = np.array(xlr_)
            nn = np.array(n_)
            idx_finite = np.isfinite(nn)
            f_finite = interpolate.interp1d(xx[idx_finite], nn[idx_finite])
            nn = f_finite(xx)
            n_ = nn.tolist()

            # solutions should not be greater than one or less than zero
            for i in range(len(n_)):
                if n_[i] < 0.:
                    n_[i] = np.round(0.0, 2)
                if n_[i] > 0.9999:
                    n_[i] = np.round(1.0, 2)

if __name__ == "__main__":
    # qsq2, c^2, g, filename
    t1 = time.time()
    master(0.165, 6.35, 1.135, 'results_hy-0.2.csv')
    t2 = time.time()

    hours = (t2 - t1)/3600
    print(str(hours) + ' hours')
