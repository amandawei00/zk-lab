# integration over zeros of Bessel function for oscillatory integrands
import numpy as np
from scipy import special
from scipy.integrate import quad

j0_zeros = special.jn_zeros(0, 1001)
def osc_intg(f, tol=1000, *args):
    x = args[0]
    k = args[1]
    sum_ = quad(f, 0, j0_zeros[0]*k, epsabs=0.0, epsrel=0.05, args=(x, k))[0]

    for i in range(tol):
        sum_ += quad(f, j0_zeros[i] * k, j0_zeros[i+1] * k, epsabs=0.0, epsrel=1.e-4, args=(x, k))[0]
    return sum_
