# tests various implementations of fourier transform of 1 - N(r,Y)
import sys
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.fft import fft, ifft, fftfreq, fftshift
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

sys.path.append('python-scripts')
sys.path.append('../bk/')
from bk_interpolate import N

y = 0
x = 0.01
i = 1j

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
bk = N('../bk/results/bk_MV1.csv')
n  = 0
m  = 0
num = 2000

def f1(r, y):
    return 1 - bk.n(r, y)

def q(t):
    frac = gamma(0.5 * (n - m - i * t + 1))/gamma(0.5 * (n + m + i * t + 1))
    return (1/2/np.pi) * np.power(2, -m - i * t) * frac

def fft1(rho):
    return (1/2/np.pi) * np.exp((1 + m) * rho) * f1(np.exp(rho), y)

p     = [fft1(k) for k in np.linspace(-100,100,num)]
pfft  = fft(p, n=num)
pfftx = fftfreq(num)

pfft  = fftshift(pfft)
pfftx = fftshift(pfftx)
phi = interp1d(pfftx, pfft)

def f2(t):
    return q(t) * phi(t)
print(pfftx)
g_    = [f2(k) for k in np.linspace(-0.49, 0.49, num)]
gifft = ifft(g_, n=num)
gifftx = np.arange(0, num, 1) # n/Fs, 1/Fs, Fs: sampling freq, 2x max freq in signal

# gifft = fftshift(gifft)
# gifftx = fftshift(gifftx)
print(gifftx)
gg = interp1d(gifftx, gifft)
def g(k):
    kap = np.log(k)
    return 2 * np.pi * np.exp((m - 1) * kap) * gg(kap) 

k_range = np.logspace(-1, 1, 100)
foo     = [g(k_range[i]) for i in range(len(k_range))]

plt.plot(k_range, foo)
plt.show() 
