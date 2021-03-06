import numpy as np
from scipy.interpolate import interp2d
from scipy.integrate import quad
from scipy.special import kn
import pandas as pd
import csv
import matplotlib.pyplot as plt

# Units in GeV
alpha = 1./137   # FIND VALUE (EM coupling)

light = 3  # q flavors CHECK 
ml    = 0.14 # q masses in GeV
el    = [2/3, -1/3, -1/3] # q charge

bk       = ''

# q20, gamma, ec = 0., 0., 0.
# bk_ is interpolated object
def mv(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * r * r * q20
    t2 = 1/(0.241 * r) + np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))

def mvg(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * np.power(r * r * q20, gamma)
    t2 = 1/(0.241 * r) + np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))

def mve(r, x):
    q20 = np.power(x0/x, lamb)
    t1 = 0.25 * r * r * q20
    t2 = 1/(0.241 * r) + ec * np.exp(1)
    return 1 - np.exp(-t1 * np.log(t2))

# (transverse) wave function for splitting of photon to quark-antiquark dipole
def psi_t2(z, r, *args): # args = [qsq2]
    coeff = (6 * alpha)/(4 * np.pi * np.pi)
    s     = 0

    for i in range(light):   # summing over light flavors
        eta2 = eta_squared(z, ml, args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(kn(0, eta * r), 2)
        k12  = np.power(kn(1, eta * r), 2)

        t1   = (z * z + (1 - z) * (1 - z)) * eta2 * k12
        t2   = ml * ml * k02
        s    += el[i] * el[i] * (t1 + t2)
    return coeff * s

# (longitudinal) wave function for splitting of photon to quark-antiquark dipole
def psi_l2(z, r, *args): # args = [qsq2]

    coeff = (6 * alpha) / (4 * np.pi * np.pi)
    s     = 0

    for i in range(light):
        eta2 = eta_squared(z, ml, args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(kn(0, eta * r), 2)
        s    += 4 * args[0] * np.power(z * (1 - z), 2) * k02 * el[i] * el[i]
    return coeff * s
        
def eta_squared(z, m_f, qsq2):
    return z * (1 - z) * qsq2 + m_f * m_f

def t_integral(z, *args): # args = [qsq2, x]
    # x = args[1] * (1 + 4 * ml * ml/args[0])
    x = args[1]
    if bk == 'mv':
        m = lambda r_: r_ * psi_t2(z, r_, args[0]) * mv_f(r_, x)
    elif bk == 'mvg':
        m = lambda r_: r_ * psi_t2(z, r_, args[0]) * mvg_f(r_, x)
    elif bk == 'mve':
        m = lambda r_: r_ * psi_t2(z, r_, args[0]) * mve_f(r_, x)

    return 2 * np.pi * quad(m, 1e-6, 1e2, epsabs=1e-5, epsrel=0.0)[0]

def l_integral(z, *args): # *args = [qsq2, x]
    # x = args[1] * (1 + 4 * ml * ml/args[0])
    x = args[1]
    if bk == 'mv':
        m = lambda r_: r_ * psi_l2(z, r_, args[0]) * mv_f(r_, x)
    elif bk == 'mvg':
        m = lambda r_: r_ * psi_l2(z, r_, args[0]) * mvg_f(r_, x)
    elif bk == 'mve':
        m = lambda r_: r_ * psi_l2(z, r_, args[0]) * mve_f(r_, x)
    return 2 * np.pi * quad(m, 1e-6, 1e2, epsabs=1e-5, epsrel=0.0)[0]

def t_xsection(x, qsq2):
    return 2 * 2.5681 * sigma * quad(t_integral, 0., 0.5, epsabs=1e-4, epsrel=0.0, args=(qsq2, x))[0]
                
def l_xsection(x, qsq2):
    return 2 * 2.5681 * sigma * quad(l_integral, 0., 0.5, epsabs=1e-4, epsrel=0.0, args=(qsq2, x))[0]

def fl(x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * l_xsection(x, qsq2)
 
def f2(x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2) + l_xsection(x, qsq2))

def reduced_x(x, qsq2, cme):
    global q20
    
    y = qsq2/(cme * cme * x)
    d = 1 + (1 - y) * (1 - y)

    a = f2(x, qsq2)
    b = (y * y/d) * fl(x, qsq2)
    c = a - b
    return c

def set_var(x0_, lamb_, gamma_, ec_, sigma_, bk_):
    global x0, lamb, gamma, ec, sigma, bk, mv_f, mvg_f, mve_f

    x0     = x0_
    lamb   = lamb_
    gamma  = gamma_
    ec     = ec_
    sigma  = 2 * sigma_
    bk     = bk_

    x_ = np.logspace(-5, -2, 400)
    r_ = np.logspace(-6, 2, 400)

    mv_grid  = [mv(rr, xx) for rr in r_ for xx in x_]
    mvg_grid = [mvg(rr, xx) for rr in r_ for xx in x_]
    mve_grid = [mve(rr, xx) for rr in r_ for xx in x_]

    mv_f  = interp2d(r_, x_, mv_grid, kind='cubic')
    mvg_f = interp2d(r_, x_, mvg_grid, kind='cubic')
    mve_f = interp2d(r_, x_, mve_grid, kind='cubic')


if __name__ == '__main__':

    res = 'interpolate_tester_interpolated.csv'
    set_var(6.51708e-5, 0.26569, 1.0699023, 5., 9.6079423, 'mvg')

    qsq = [1.5, 8.5, 27., 200.]
    with open(res, 'w') as f:
        writer = csv.writer(f, delimiter='\t')    
        writer.writerow(['q2', 'x', 'redx'])
        for i in qsq:
            print(i)
            for j in np.logspace(-5, -2, 20):
                writer.writerow([i, j, reduced_x(j, i, 319.)])
