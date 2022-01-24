cimport cython
import numpy as np
cimport numpy as cnp
from libc.math cimport pow

cdef double eta2(double z, double m_f, double qsq2):
    return z * (1 - z) * qsq2 + pow(m_f, 2)



