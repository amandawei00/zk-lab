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
A  = 197
nc = 3
cf = (nc * nc - 1)/(2 * nc)
s  = 1


j0_zeros = special.jn_zeros(0, 1001)
################################################################## distributions
# integrates oscillatory bessel function by summing over zeros of Bessel function
def integrate(f, tol, *args):
    k = args[0]
    # err  = 10.
    # i    = 0
    sum_ = quad(f, 0, j0_zeros[0] * k, epsabs=0.0, epsrel=0.05, args=(k,))[0]
    # while np.abs(err) > tol:
    #     err  = quad(f, j0_zeros[i], j0_zeros[i+1], epsabs=0.0, epsrel=0.05, args=(args[0],))[0]
    #     sum_ += err
    #     i    += 1
    for i in range(1000):
        sum_ += quad(f, j0_zeros[i] * k, j0_zeros[i+1] * k, epsabs=0.0, epsrel=1.e-4, args=(k,))[0]
    return sum_

def gq1_integrand(u, *args):
    k = args[0]
    a = u * u / (k * k)
    jac = a / u
    bes = special.j0(u)
    lap = 0.5 * np.exp(-0.25 * a) * (2 - 0.5 * a)
    return jac * bes * lap

# a1
def qg1(x, k_):
    prefac = nc * s_/(4 * np.pi * np.pi * np.pi)
    res = integrate(gq1_integrand, 1.e-50, k_)
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
f = [qg1(0.01, i)/5.009 for i in k]

plt.xlim(1.e-2, 1.e2)
plt.ylim(1.e-6, 1.e0)
plt.xscale('log')
plt.yscale('log')

plt.xlabel('k [GeV]')
plt.ylabel('TMD/transverse nuclear area, [GeV^-2]')
plt.title('qg(1)')
plt.plot(k, f)
plt.show()
