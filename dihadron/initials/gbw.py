from numpy import exp
import numpy as np
from sympy import symbols, diff, lambdify
import matplotlib.pyplot as plt

# qsq0 [GeV^2] = 1 [GeV^2]: initial saturation scale
# x0           = 0.01     : onset of small-x
# gamma        = 0.3      : 
def s(x, r_, x0=0.01, qsq0=1.0, gamma=0.3):
    qsq = qsq0 * np.power(x0/x, gamma)
    return exp(-(1/4) * qsq * r_ * r_)

'''
r = np.logspace(-2 ,2, 50)
y = [s(0.01, r[i]) for i in range(len(r))]

plt.xscale('log')
plt.yscale('log')
plt.plot(r,y)
plt.show()
'''
