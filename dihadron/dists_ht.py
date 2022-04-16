import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.special import j0
from scipy.fft import fht
from scipy.signal import convolve
from integrator import osc_intg

import sys
sys.path.append('./initials/')


# init 
initial = 'mv'
x0      = 0.01
g       = 0.3

if initial == 'mv':
    import mv as pig
    qsq0 = 0.6
elif initial == 'gbw':
    import gbw as pig
    qsq0 = 1.

def qsq(x):
    return qsq0 * np.power(x0/x, g)

q2      = qsq(x0)
# function parameters
n       = 399
r1      = 3.e-6
r2      = 60.e0
xr1     = np.log(r1)
xr2     = np.log(r2)
hr      = (xr2 - xr1) / n

xlr_    = [xr1 + i * hr for i in range(n + 1)]
r_      = np.exp(xlr_)

klr_    = [1/xlr_[i] for i in range(len(xlr_))]
k_      = [1/r_[i] for i in range(len(r_))]

# parameters
A  = 197
nc = 3
cf = (nc * nc - 1)/(2 * nc)
s_ = 1
c0 = cf * s_/(4 * np.pi * np.pi * np.pi)
c1 = nc * s_/(4 * np.pi * np.pi * np.pi)
# args = [x, k]

def h(f, g, c0):
    # hankel transform functions f and g
    f_ = fht(f, dln=hr, mu=0) # dln=np.exp(hr)?
    g_ = fht(g, dln=hr, mu=0)

    # convolve f and g
    c  = convolve(f_, g_, mode='same')
    return -(1/(2 * np.pi)) * interp1d(k_, c)

def qg1_f(x, r):
    return 1 - pig.s(x, r)

def qg1_g(x, r):
    return c1

def qg2_f(x, r):
    return - np.log(pig.s(x, r))

def qg2_g(x, r):
    m = pig.s(x, r)
    return ((1 - np.power(m, nc/cf)) * m)/(-np.log(m))

def gg1_f(x, r):
    return 1 - pig.s(x, r)

def gg1_g(x, r):
    return c1 * pig.s(x, r)

def adj_f(x, r):
    return 1 - np.power(pig.s(x, r), nc/cf)

def adj_g(x, r):
    return c1

def ww_f(x, r):
    return -np.log(pig.s(x, r))

def ww_g(x, r):
    m = pig.s(x, r)
    return c0 * (1 - np.power(m, nc/cf))/(-m)

def dist(f, x, k_):
    res = osc_intg(f, 1000, x, k_)
    print('F(k  = ' + str(k_) + ') = ' + str(res))
    return res

k = np.logspace(-2, 2, 100)
qg1_ = [dist(qg1, 0.01, i) for i in k]
qg2_ = [dist(qg2, 0.01, i) for i in k]
gg1_ = [dist(gg1, 0.01, i) for i in k]
adj_ = [dist(adj, 0.01, i) for i in k]
ww_  = [dist(ww,  0.01, i) for i in k]
gg6_ = [dist(gg6, 0.01, i) for i in k]

plt.xlim(1.e-2, 1.e2)
plt.ylim(1.e-6, 1.e0)
plt.xscale('log')
plt.yscale('log')

lw = 1
plt.xlabel('k [GeV]')
plt.ylabel('TMD/transverse nuclear area, [GeV^-2]')
plt.title('small-x TMDs at Y=0 (x=0.01)')
plt.plot(k, qg1_, 'black', linestyle='-', linewidth=lw, label='qg(1)')
plt.plot(k, qg2_, 'blue', linestyle='--', linewidth=lw, label='qg(2)')
plt.plot(k, gg1_, 'orange', linestyle='-.', linewidth=lw, label='gg(1)')
plt.plot(k, adj_, 'magenta', linestyle=':', linewidth=lw, label='adj')
plt.plot(k, gg6_, 'green', linestyle='--', linewidth=lw, label='gg(6)')
plt.plot(k, ww_, 'cyan', linestyle='-', linewidth=lw, label='WW')
plt.legend()
plt.show()

