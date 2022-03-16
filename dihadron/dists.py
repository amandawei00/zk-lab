import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import special
from integrator import osc_intg

import sys
sys.path.append('./initials/')

# init
initial = 'mv'
x0      = 0.01
g       = 0.3 

if initial is 'mv':
    import mv as pig
    qsq0 = 0.6
elif initial is 'gbw':
    import gbw as pig
    qsq0 = 1.

def qsq(x):
    return qsq0 * np.power(x0/x, g)

q2 = qsq(x0)

# parameters
A  = 197
nc = 3
cf = (nc * nc - 1)/(2 * nc)
s_ = 1

# args = [x, k]
def gq1(u, *args):
    x   = args[0]
    k   = args[1]
    pre = nc * s_/(4 * np.pi * np.pi * np.pi)
    jac = u / (k * k)
    bes = special.j0(u)
    return pre * jac * bes * pig.lap(x, u/k)

def qg2(u, *args):
    x   = args[0]
    k   = args[1]
    pre = cf * s_/(4 * np.pi * np.pi * np.pi)

    jac = u / (k * k)
    bes = special.j0(u)
    mid = pig.kap(x, u/k) * (1 - np.power(pig.s(x, u/k), nc/cf))
    return pre * jac * bes * mid * pig.s(x, u/k)

def gg1(u, *args):
    x   = args[0]
    k   = args[1]
    pre = nc * s_/(4 * np.pi * np.pi * np.pi)
    jac = u / (k * k)

    bes = special.j0(u)
    return pre * jac * bes * pig.s(x, u/k) * pig.lap(x, u/k)

def adj(u, *args):
    x   = args[0]
    k   = args[1]
    pre = cf * s_/(4 * np.pi * np.pi * np.pi)
    jac = u / (k * k)

    bes = special.j0(u)
    mid = pig.lap2(x, u/k, nc/cf)
    return pre * jac * bes * mid

def ww(u, *args):
    x   = args[0]
    k   = args[1]
    pre = cf * s_/(4 * np.pi * np.pi * np.pi)
    jac = u / (k * k)

    bes = special.j0(u)
    mid = pig.kap(x, u/k) * (1 - np.power(pig.s(x, u/k), nc/cf))
    return pre * jac * bes * mid

def gg6(u, *args):
    x   = args[0]
    k   = args[1]
    pre = cf * s_/(4 * np.pi * np.pi * np.pi)
    jac = u / (k * k)

    bes = special.j0(u)
    mid = pig.kap(x, u/k) * (1 - np.power(pig.s(x, u/k), nc/cf))
    end = np.power(pig.s(x, u/k), nc/cf)

    return pre * jac * bes * mid * end    

def dist(f, x, k_):
    res = osc_intg(f, 1000, x, k_)
    print('F(k  = ' + str(k_) + ') = ' + str(res))
    return res

k = np.logspace(-2, 2, 50)
f = [dist(qg2, 0.01, i) for i in k] # /5.009

plt.xlim(1.e-2, 1.e2)
plt.ylim(1.e-6, 1.e0)
plt.xscale('log')
plt.yscale('log')

plt.xlabel('k [GeV]')
plt.ylabel('TMD/transverse nuclear area, [GeV^-2]')
plt.title('qg(1)')
plt.plot(k, f)
plt.show()
