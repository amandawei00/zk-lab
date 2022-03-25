import numpy as np

y    = 0
x0   = 0.01
qsq0 = 1.   # GeV^2
a    = 0.3  # little lambda

def qsq(x):
    return qsq0 * np.power(x0/x, a)

def s(x, r):
    q2 = qsq(x)
    return np.exp(-0.25 * q2 * r * r)

# derivative of N(x,r) = 1 - S(x,r) with respect to r
def dn(x, r):
    q2 = qsq(x)
    return 0.5 * q2 * r * np.exp(-0.25 * q2 * r * r)

# second derivative of N(x,r) with respect to r
def d2n(x, r):
    q2 = qsq(x)
    r2 = r * r

    ex = np.exp(-0.25 * q2 * r2)
    return ex * (-0.25 * q2 * q2 * r2 + 0.5 * q2)

# laplacian of N(x, r) = 1 - S(x, r)
def lap(x, r):
    return d2n(x, r) + (1/r) * dn(x, r)

# gamma and laplacian of gamma --> kappa
def gam(x, r):
    return -np.log(s(x,r))

def dgam(x, r):
    q2 = qsq(x)
    return 0.5 * q2 * r

def d2gam(x, r):
    q2 = qsq(x)
    return 0.5 * q2

def lap1(x, r):
    return d2gam(x, r) + (1/r) * dgam(x, r)

def kap(x, r):
    return lap1(x,r)/gam(x,r)

# laplacian of 1 - S(x, r)^(nc/cf)
def df(x, r, a):
    return a * np.power(s(x, r), a-1) * dn(x, r)

def d2f(x, r, a):
    t1 = a * (a - 1) * np.power(s(x, r), a-2) * dn(x, r) * dn(x, r)
    t2 = a * np.power(s(x, r), a-1) * d2n(x, r)
    return -t1 + t2

def lap2(x, r, a):
    return d2f(x, r, a) + (1/r) * df(x, r, a)
