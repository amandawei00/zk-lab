import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import kn
import csv

def psi_t2(z, r, *args): # args = [qsq2]
    coeff = (6 * alpha)/(4 * np.pi * np.pi)
    s     = 0

    for i in range(light):   # summing over light flavors
        eta2 = eta_squared(z, ml, args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(kn(0, eta * r), 2) # modified Bessel function (kind 2, order 0)
        k12  = np.power(kn(1, eta * r), 2)  # MacDonald's Function first order

        t1   = (z * z + (1 - z) * (1 - z)) * eta2 * k12
        t2   = ml * ml * k02
        s    += el[i] * el[i] * (t1 + t2)
    return coeff * s
def eta_squared(z, m_f, qsq2):
    return z * (1 - z) * qsq2 + m_f * m_f


def t_integrand(z, *args): # args = [qsq2, y]
    return r_ * psi_t2(z, r_, args[0]) * bk.n(r_, args[1])


