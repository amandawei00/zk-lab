import numpy as np
import ctypes
from ctypes.util import find_library
from ctypes import byref, CDLL, c_int, c_double, POINTER, cast, c_float
import matplotlib.pyplot as plt
import sys
import glob

cspline = CDLL('./spline_c.so')
print("shared object imported")

n = 10
m = 50

xmin = 0.0
xmax = 7

x  = np.linspace(xmin, xmax, n)
x_grid = np.linspace(xmin, xmax, m)
y = np.sin(x)

b = np.empty(n)
c = np.empty(n)
d = np.empty(n)

cspline.spline.restype = None
cspline.spline.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_int]

cspline.ispline.restype = ctypes.c_double
cspline.ispline.argtypes = [ctypes.c_double] + cspline.spline.argtypes

cspline.spline(x, y, b, c, d, n)
y_grid = [cspline.ispline(x0, x, y, b, c, d, n) for x0 in x_grid]


plt.errorbar(x,y,fmt='o')
plt.plot(x_grid, y_grid)
plt.show()


