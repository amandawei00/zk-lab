from sympy import symbols, diff, lambdify
import sympy as sp
import numpy as np
from scipy import special
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys

sys.path.append('./initials/')
from gbw import s as s
# from mv import s as s

# variables
lamb = 0.241
qsq0 = 1.0   # GeV^2
x0   = 0.01
g    = 0.3

# parameters
nc = 3
s_ = 1

# laplacians of initial conditions
x, r = symbols('x r')
gbw_ = 1 - sp.exp(-0.25 * qsq0 * (x0/x)**g * r * r)
gbw_p1 = diff(gbw_, r)
gbw_p2 = diff(gbw_p1, r)
del_gbw = gbw_p2 + (1/r) * gbw_p1
del_s = lambdify([x, r], del_gbw)

'''x, r = symbols('x r')
mv_  = 1 - sp.exp(-0.25 * qsq0 * (x0/x)**lamb * r * r * sp.log(1/(g * r) + sp.exp(1)))
mv_p1 = diff(mv_, r)
mv_p2 = diff(mv_p1, r)
del_mv = mv_p2 - (1/r) * mv_p1
del_s = lambdify([x,r], del_mv)'''

gm  = -sp.log(gbw_)
gm1 = diff(gm, r)
gm2 = diff(gm1, r)
del_gm = lambdify([x,r], gm2 + (1/r) * gm1)

# useful functions
def k(x, r_):
    return del_gm(x, r_)/gm(x, r_)
    
def gm(x, r_):
    return -np.log(s(x, r_))

j0_zeros = special.jn_zeros(0, 100)
print(j0_zeros)
################################################################## distributions
# integrates oscillatory bessel function by summing over zeros of Bessel function
def integrate(f, tol):
    err  = 10.
    sum_ = 0.
    i    = 1
    while err > tol:
        err  = quad(f, j0_zeros[i], j0_zeros[i+1], epsabs=0.0, epsrel=1.e-8)[0]
        sum_ += err
        i    += 1
    return sum_

# a1
def qg1(x, k_):
    prefac = nc * s_/(2 * 2 *  np.pi * np.pi * np.pi)
    f = lambda u: u * special.j0(k_ * u) * del_s(x, u) 
    res = integrate(f, 1.e-10)
    print('F(k  = ' + str(k_) + ') = ' + str(res))
    return res

# a2
# def qg2(x, k_):
# a3
# def gg1(x, k_):
# a4
# def adj(x, k_):
# a5
# def ww(x, k_):
# a6
# def gg6(x, k_):

k = np.logspace(-2, 2, 50)
f = [qg1(0.01, i) for i in k]

plt.xlim(1.e-2, 1.e2)
plt.ylim(1.e-8, 1.e0)
plt.xscale('log')
plt.yscale('log')
plt.plot(k, f)
plt.show()
