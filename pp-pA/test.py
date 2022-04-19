# tests various implementations of fourier transform of 1 - N(r,Y)
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

sys.path.append('python_scripts')
from bk_interpolate import N
import lhapdf as pdf

# direct, straightforward implementation
'''n1  = N('../bk/results/bk_MV1.csv')
n2 = N('../bk/results/bk_MVe1.csv')
n3 = N('../bk/results/bk_MVg1.csv')

print('data loaded...')
y = 0
x = 0.01

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
