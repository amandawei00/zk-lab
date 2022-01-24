from sympy import symbols, diff, lambdify
import sympy as sp
import numpy as np
from scipy.special import j0
import sys

sys.path.append('./initials/')
from gbw import s as s

# variables
lamb = 0.241
qsq0 = 1.0   # GeV^2
x0   = 0.01
g    = 0.3

# parameters
nc = 3
s_ = 1

# laplacians of initial conditions
x, r = symbols('x r')
gbw_ = sp.exp(-0.25 * qsq0 * (x0/x)**lamb * r * r)
gbw_p1 = diff(gbw_, r)
gbw_p2 = diff(gbw_p1, r)
del_gbw = gbw_p2 - (1/r) * gbw_p1
del_s = lambdify([x, r], del_gbw)

'''x, r = symbols('x r')
mv_  = 1 - sp.exp(-0.25 * qsq0 * (x0/x)**lamb * r * r * sp.log(1/(g * r) + sp.exp(1)))
mv_p1 = diff(mv_, r)
mv_p2 = diff(mv_p1, r)
del_mv = mv_p2 - (1/r) * mv_p1
del_s = lambdify([x,r], del_mv)'''

gm = -sp.log(gbw_)
gm = diff(gm, r)
gm = diff(gm1, r)
del_gm = lambdify([x,r], gm2 - (1/r) * gm1)
################################################################### useful functions
def k(x, r_):
    return del_gm(x, r_) - gm(x, r_)
    
def gm(x, r_):
    return -np.log(s(x, r_))
################################################################## distributions

# a1
def qg1(x, k_):
    prefac = nc * s_/(2 * np.pi * np.pi)
    integrand = 
# a2
def qg2(x, k_):
# a3
def gg1(x, k_):
# a4
def adj(x, k_):
# a5
def ww(x, k_):
# a6
def gg6(x, k_):

