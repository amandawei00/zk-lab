from scipy.signal import convolve
from scipy.special import gamma
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

# testing Hankel transform on function of the form (r^n)exp^{-a^2 r^2}
# HT of (r^n)exp{-a^2 r^2} from table:
def f(s, n, a):
    ex = np.exp(-s * s/(4 * a * a))
    return np.power(s/(2 * a * a), n) * ex/(2 * a * a)

# HT of exp{- epsilon * r} * (r^n) * exp{-a^2 r^2}
def g(n, a, ep):
    ss  = np.linspace(0, 100, 300)
    num = lambda x: 2 * ep * np.power(2 * x, n) * gamma(n + 3/2)
    den = lambda x: np.power(x * x + ep * ep, n + 3/2) * np.sqrt(np.pi)
    ex  = lambda x: (1/2/a/a) * np.exp(-x * x/(4 * a * a))
    g1  = [num(ss[i])/den(ss[i]) for i in range(len(ss))]
    g2  = [ex(ss[i]) for i in range(len(ss))]
    g3  = convolve(g1, g2, mode='same')
    return interp1d(ss, g3, kind='cubic')

ep = 0.1
n  = 0    # order of Hankel Transform
a  = 0.5  # some constant value
s  = np.linspace(0, 100, 300)
ff = [f(s[i], n, a) for i in range(len(s))]
g_ = g(n, a, ep)
gg = [g_(s[i]) for i in range(len(s))]

plt.plot(s, ff, label='HT table solution')
plt.plot(s, gg, label='solution from regularization')
plt.xscale('log')
plt.legend()
plt.show()
