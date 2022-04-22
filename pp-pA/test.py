# tests various implementations of fourier transform of 1 - N(r,Y)
import sys
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.integrate import quad
from scipy.fft import fft, ifft, fftfreq, fftshift
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

sys.path.append('python-scripts')
sys.path.append('../bk/')
from bk_interpolate import N

#y = 0
#x = 0.01
#i = 1j

# direct, straightforward implementation
'''
n1  = N('../bk/results/bk_MV1.csv')
n2 = N('../bk/results/bk_MVe1.csv')
n3 = N('../bk/results/bk_MVg1.csv')

print('data loaded...')
k_range = np.logspace(-1, 1, 100)

# fundamental representation
t1 = time.time()
udgf1 = [n1.udg_f(x, k) for k in k_range]
t2 = time.time()
udgf2 = [n2.udg_f(x, k) for k in k_range]
t3 = time.time()
udgf3 = [n3.udg_f(x, k) for k in k_range]
t4 = time.time()

print('time taken: ' + str((t2-t1)/60) + ' min, ' + str((t3-t2)/60) + ' min, ' + str((t4-t3)/60) + ' min')

plt.xlim(1.e-1, 1.e1)
plt.ylim(1.e-4, 1.e2)
plt.xscale('log')
plt.yscale('log')

plt.xlabel('k')
plt.ylabel('Fourier Transform of S(r, Y=0) = 1 - N(r, Y=0)')
plt.plot(k_range, udgf1, '-', color='black', label='MV initial conditions')
plt.plot(k_range, udgf2, '--', color='blue', label='MVe initial conditions')
plt.plot(k_range, udgf3, ':', color='magenta', label='MVg initial conditions')
plt.legend()
plt.show()
'''


# implementation with Hankel Transforms
t1 = time.time()
bk = N('../bk/results/bk_MV1.csv')
print(str(time.time()-t1))
n  = 0
m  = 0

from scipy.special import j0

def f(k, y):
    integrand = lambda r: r * j0(k*r) * (1 - bk.n(r,y))
    return quad(integrand, 3.e-6, 60.)[0]

# rrange = np.linspace(0, 10, 100)
# bk_ = [bk.n(rrange[i], 0) for i in range(len(rrange))]
# plt.plot(rrange, bk_)
# plt.show()

krange = np.linspace(0.1, 10, 25)
# y = [f(b) for b in krange]
y = []
y1 = []
for i in range(len(krange)):
    t1 = time.time()
    y.append(f(krange[i], 0))
    y1.append(f(krange[i], 0.1))
    print('integral took ' + str((time.time()-t1)) + ' s')

# y = [f(krange[i]) for i in range(len(krange))]
print(krange)
print(y)
print(y1)
plt.yscale('log')
plt.plot(krange, y)
plt.plot(krange, y1)
plt.show()

# r_ = np.linspace(0, 10, 25)
# y = [integrand(i) for i in r_]


def f1(r, y):
    return 1 - bk.n(r, y)

def q(t):
    frac = gamma(0.5 * (n - m - i * t + 1))/gamma(0.5 * (n + m + i * t + 1))
    return (1/2/np.pi) * np.power(2, -m - i * t) * frac

def fft1(rho):
    return (1/2/np.pi) * np.exp((1 + m) * rho) * f1(np.exp(rho), y)

'''
num   = 400
T     = 1./8
x1    = np.linspace(0.0, num * T, num)
y1    = [fft1(x_) for x_ in x1]

x1f   = np.linspace(0.0, 1.//(2.0 * T), num//2)
phi   = interp1d(x1f, y1f)

def f2(t):
    return q(t) * phi(t)

x2 = np.linspace(0, 3.5, 100)
y2 = [f2(x2[i]) for i in range(len(x2))]

y2f = ifft(y2, n=num)
x2f = np.arange(0, num/40, 1/40)
gg  = interp1d(x2f, y2f)

def g(k):
    kap = np.log(k)
    return 2 * np.pi * np.exp((m - 1) * kap) * gg(kap)

x_range = np.logspace(-1, 1, 100)
y_range = [g(k) for k in x_range]

plt.plot(x_range, y_range)
plt.show()

# ht version without use of fft just to check

def phi(t):
    f = lambda rho: np.sin(t * rho) * np.exp((1 + m) * rho) * f1(np.exp(rho), 0)
    return (1/2/np.pi) * quad(f, -50, 50)[0]

def goo(kap):
    f = lambda t: np.sin(-kap * t) * q(t) * phi(t)
    return 2 * np.pi * np.exp((m - 1) * kap) * quad(f, -50, 50)[0]

krange = np.logspace(-1, 1, 100)
y = [goo(np.exp(kap)) for kap in krange]

plt.plot(krange, y)
plt.show()'''
