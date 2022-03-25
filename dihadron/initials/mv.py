import numpy as np

y    = 0
x0   = 0.01
qsq0 = 0.6    # GeV^2
a    = 0.3   # little lambda
lamb = 0.241 # Lambda QCD

def qsq(x):
    return qsq0 * np.power(x0/x, a)

def s(x, r):
    q2 = qsq(x)
    return np.exp(-0.25 * q2 * r * r * np.log(1/(lamb * r) + np.exp(1)))

# derivative of N(x,r) = 1 - S(x,r) with respect to r
def dn(x, r):
    q2   = qsq(x)
    beta = np.exp(1) + 1/(lamb * r)

    t1   = -np.power(beta, -0.25 * q2 * r * r)
    t2   = 0.25 * q2/(lamb * beta)
    t3   = 0.5 * q2 * r * np.log(beta)
    return t1 * (t2 - t3)

# second derivative of N(x,r) with respect to r
def d2n(x, r):
    q2   = qsq(x)
    beta = np.exp(1) + 1/(lamb * r)
    gbr  = lamb * beta * r

    t1   = np.power(beta, -0.25 * q2 * r * r)
    t2   = 0.25 * q2/(gbr * gbr)
    t3   = 0.50 * q2/(gbr)
    t4   = 0.50 * q2 * np.log(beta)
    t5   = r * (t3/2 - t4)

    return -t1 * (t2 + t3 - t4) - t1 * (t5 * t5)

# laplacian of N(x, r) = 1 - S(x, r)
def lap(x, r):
    return d2n(x, r) + (1/r) * dn(x, r)

# gamma and laplacian of gamma --> kappa
def gam(x, r):
    return -np.log(s(x, r))

def dgam(x, r):
    q2 = qsq(x)
    beta = np.exp(1) + 1/(lamb * r)
    db = -1/(lamb * r * r)

    t1 = 0.50 * q2 * r * np.log(beta)
    t2 = 0.25 * q2 * r * r * (1/beta) * db
    return t1 + t2

def d2gam(x, r):
    q2 = qsq(x)
    beta = np.exp(1) + 1/(lamb * r)
    db = -1/(lamb * r * r)
    d2b = 2/(lamb * r * r * r)

    t1 = 0.5 * q2 * np.log(beta)
    t2 = q2 * r * (1/beta) * db
    t3 = 0.25 * q2 * r * r * ((1/beta) * d2b - (1/(beta * beta)) * db)

    return t1 + t2 + t3

def lap1(x, r):
    return d2gam(x, r) + (1/r) * dgam(x,r)

def kap(x, r):
    return lap1(x,r)/gam(x,r)

# laplacian of 1 - S(x, r)^(nc/cf)
def df(x, r, a):
    return a * np.power(s(x, r), a-1) * dn(x, r)

def d2f(x, r, a):
    t1 = a * (a - 1) * np.power(s(x, r), a-2) * dn(x,r) * dn(x,r)
    t2 = a * np.power(s(x, r), a-1) * d2n(x,r)
    return -t1 + t2

def lap2(x, r, a):
    return d2f(x, r, a) + (1/r) * df(x, r, a)
