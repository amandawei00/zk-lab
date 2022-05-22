from scipy import LowLevelCallable as llc
import scipy.integrate as intg
from numba.types import intc, voidptr, CPointer, float64
from numba import cfunc
import t
import subprocess
import numpy as np
import ctypes
from ctypes import c_void_p, c_int, POINTER, c_uint64

@cfunc(float64(intc, CPointer(float64), voidptr))
def f1_wrap(n, x, udata):
    return t.f1(x[0], udata[0])

udata = np.array([2])
# f1 = llc.from_cython(t, 'f1', user_data=CPointer(float64))
# print(intg.quad(f1, 0., 4., args=(2,)))
print(intg.quad(LowLevelCallable(f1_wrap, (2,)), 0., 4.))

#------------------------------------------------------------------------------------
# subprocess.run(['python3', 'setup.py', 'build_ext', '--inplac'])

# def sample(x, y):
#     return x * x + 2 * y

# single_intg = llc.from_cython(t, 'py_f3')
# print(intg.quad(single_intg, 0, 4))

# print("python integral")
# print(intg.nquad(sample, [[0, 4], [0, 2]]))

# double_intg = llc.from_cython(t, 'py_f2', signature='double (int, double *)')
# print("cython integral")

# x * x + 2 * y, xx=[x,y]. assuming x is first
# print(intg.dblquad(double_intg, 0, 6, 0, 2, epsabs=0.0, epsrel=0.5))
# print(intg.dblquad(double_intg, 0, 2, 0, 6, epsabs=0.0, epsrel=0.5))
#print(intg.nquad(double_intg, [[0, 4],[0,2]]))

# print(t.py_f(np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.intc), 7))
# t.set_abc(3,2, 5)
# t.print_abc()
# t.set_abc(2, 2, 2)
# t.print_abc()
