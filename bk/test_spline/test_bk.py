import numpy as np
import ctypes
from ctypes.util import find_library
from ctypes import byref, CDLL, c_int, c_double, POINTER, cast, c_float
import matplotlib.pyplot as plt
import sys
import glob
from numpy import sqrt, exp
import pandas as pd
# variables
n = 399     # number of r points to be evaluated at each evolution step in Y
r1 = 1.e-6  # limits of r
r2 = 1.e2

xr1 = np.log(r1)
xr2 = np.log(r2)

hr = (xr2 - xr1) / n

hy = 0.2
ymax = 16.
y = np.arange(0.0, ymax, hy)

# Arrays for N and r in N(r), evaluated at some rapidity Y (including next step N(r,Y) in the evolution
xlr_ = [xr1 + i * hr for i in range(n + 1)]
r_ = np.exp(xlr_)
n_ = []

# parameters
nc = 3        # number of colors
nf = 3        # number of active flavors
lamb = 0.241  # lambda QCD (default)

beta = (11 * nc - 2. * nf)/(12 * np.pi)
afr = 0.7     # frozen coupling constant (default)

c2, gamma, qs02, ec = 0. , 0., 0., 0.   # fitting parameters
e  = np.exp(1)

def mv(r):
    xlog = np.log(1/(lamb * r) + ec * np.exp(1))
    xexp = np.power(qs02 * r * r, gamma) * xlog/4.0
    return 1 - np.exp(-xexp)

cspline = CDLL('./spline_c.so')
print("shared object imported")

cspline.spline.restype = None
cspline.spline.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_int]

cspline.ispline.restype = ctypes.c_double
cspline.ispline.argtypes = [ctypes.c_double] + cspline.spline.argtypes

def nfunc(qlr, xlr_, n_, coeff1, coeff2, coeff3):
    x = 0.0
    if qlr < xr1:
        x = n_[0] * exp(2 * qlr)/(1.e-6 * 1.e-6)
    elif qlr >= xr2:
        x = 1.
    else:
        x = cspline.ispline(qlr, xlr_, n_, coeff1, coeff2, coeff3, n)

    if x < 0.: return 0.0
    if x > 1.: return 1.0

    return x

def alphaS(rsq):
    if rsq > rfr2:
        return afr
    else:
        xlog = log((4 * c2)/(rsq * lamb * lamb))
        return 1/(beta * xlog)

def find_r1(r, z, thet):
    r12 = (0.25 * r * r) + (z * z) - (r * z * cos(thet))
    return sqrt(r12)


def find_r2(r, z, thet):
    r22 = (0.25 * r * r) + (z * z) + (r * z * cos(thet))
    return sqrt(r22)


# kernel
def k(r, r1_, r2_):
    if (r1_ < 1e-20) or (r2_ < 1e-20):
        return 0
    else:
        rr = r * r
        r12 = r1_ * r1_
        r22 = r2_ * r2_

        t1 = rr/(r12 * r22)
        t2 = (1/r12) * (alphaS(r12)/alphaS(r22) - 1)
        t3 = (1/r22) * (alphaS(r22)/alphaS(r12) - 1)

        prefac = (nc * alphaS(rr))/(2 * M_PI * M_PI)
        return prefac * (t1 + t2 + t3)


def f_kernel(theta, xlr):  # xx = [theta, a]
    z = exp(xlr)
    r1_ = find_r1(r0, z, theta)
    r2_ = find_r2(r0, z, theta)

    return z * z * k(r0, r1_, r2_)

def f_split(theta, xlr):

    z = exp(xlr)
    r1_ = find_r1(r0, z, theta)
    r2_ = find_r2(r0, z, theta)

    nr1 = nfunc(log(r1_))
    nr2 = nfunc(log(r2_))

    return z * z * k(r0, r1_, r2_) * (nr1 + nr2)

# combined integrand
def f_combined(theta, xlr):

    z = exp(xlr)
    r1_ = find_r1(r0, z, theta)
    r2_ = find_r2(r0, z, theta)

    xlr1 = log(r1_)
    xlr2 = log(r2_)

    nr1 = nfunc(xlr1)
    nr2 = nfunc(xlr2)

    return z * z * k(r0, r1_, r2_) * (nr1 + nr2 - n0 - nr1 * nr2)

prev = 0.5  # previous step in rapidity
f    = '../results/bk_MV1.csv'

# read from file
df     = pd.read_csv(f, delimiter='\t', header=0, index_col=None)
r_prev = df.loc[df['y']==prev][['r']]
n_prev = df.loc[df['y']==prev][['N(r,Y)']]

b = np.empty(n)
c = np.empty(n)
d = np.empty(n)
print(r_prev)
print(n_prev)
cspline.spline(r_prev.to_numpy(), n_prev.to_numpy(), b, c, d, n)

x      = np.logspace(xr1, xr2, n)
x_grid = np.logspace(xr1, xr2, 600)
y      = np.array([nfunc(xlr, r_prev, n_prev, b, c, d) for xlr in x])
y_grid = [nfunc(xlr, r_prev, n_prev, b, c, d) for xlr in x_grid]

# ratio = [y_grid[i]/ for i in range(len(y_grid))]
plt.xscale('log')
plt.errorbar(x,y,fmt='o', markersize=1.)
plt.plot(x_grid, y_grid)
# print(ratio)
# plt.plot(x_grid, ratio)
plt.show()


