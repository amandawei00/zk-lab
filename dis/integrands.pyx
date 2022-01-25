cimport cython
cimport scipy.special.cython_special as cspec
import numpy as np
cimport numpy as cnp
from libc.math cimport pow, M_PI, sqrt

cdef double alpha =1./137.
cdef double sig   = 1.
cdef double x0    = 0.01
cdef double qsq2  = 0.4
cdef double lamb = 0.241

flavors = [1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6]
mf = [0.002, 0.0045, 1.270, 0.101, 172, 5., 0.002, 0.0045, 1.270, 0.101, 172., 5.]
ef = [2/3, -1/3, 2/3, -1/3, 2/3, -1/3, -2/3, 1/3, -2/3, 1/3, -2/3, 1/3]

cpdef void set_var(qsq2_):
    global qsq2
    qsq2 = qsq2_

cdef double eta2(double z, double m_f, double qsq2):
    return z * (1 - z) * qsq2 + pow(m_f, 2)

cdef double psi_t2(int n, double *xx): # xx = [z, r]
    cdef double coeff = (6 * alpha)/(4 * M_PI * M_PI)
    cdef double summ = 0
    for i in range(len(flavors)):
        e2 = eta2(xx[0], mf[i], qsq2)
        eta = sqrt(e2)
        k02 = pow(cspec.kn(0, eta * xx[1]), 2)
        k12 = pow(cspec.kn(1, eta * xx[1]), 2)
        t1 = (pow(xx[0], 2) + pow(1 - xx[0], 2)) * e2 * k12
        t2 = pow(mf[i], 2) * k02
        summ += pow(ef[i], 2) * (t1 + t2)
    return coeff * summ

cdef double psi_l2(int n, double *xx): # xx=[z, r]
    cdef double codeff = (6 * alpha) / (4 * M_PI * M_PI)
    cdef double summ = 0.

    for i in range(len(flavors)):
        e2 = eta2(xx[0], mf[i], qsq2)
        eta = sqrt(e2)
        k02 = pow(cspec.kn(0, eta * xx[1]), 2)
        summ += 4 * qsq2 * pow(xx[0] * (1 - xx[0]), 2) * k02 * pow(ef[i], 2)

    return coeff * summ


