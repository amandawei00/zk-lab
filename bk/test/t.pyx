cimport cython
from cython.view cimport array as cvarray
from cpython.pycapsule cimport (PyCapsule_New, PyCapsule_GetPointer)
import numpy as np
cimport numpy as cnp
from libc.math cimport log
from libc cimport int
cnp.import_array()

cdef int a = 0
cdef int b = 0
cdef int c = 0
cdef double d = log(10)

cdef extern from "clyde.c":
    int func(int *x, int n)
    double func2(int n, double *xx)
    double func3(double x)

# def py_f(np.ndarray[int, ndim=1, mode='c'] x, int n):
cpdef py_f(cnp.ndarray[int, ndim=1, mode='c'] x, int n):
    # cdef int [:] x_ = x
    # cdef int[:] x_ = np.ascontiguousarray(x, dtype=int)
    return func(<int*> x.data, n)

cdef double py_f2(int n, double *xx):
    return func2(n, xx)

cdef double py_f3(double x):
    return func3(x)

cpdef void set_abc(int ant, int baboon, int clown):
    global a, b, c

    a = ant
    b = baboon
    c = clown

cpdef void print_abc():
    print("a = " + str(a))
    print("b = " + str(b))
    print("c = " + str(c))
