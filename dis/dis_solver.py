import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import kn
import pandas as pd
import csv

from os.path import exists

sys.path.append('../bk/')
from bk_interpolate import N

# Units in GeV
alpha = 1./137   # FIND VALUE (EM coupling)
x0    = 0.01
lamb  = 0.241  # lambda_QCD (GeV)

light = 3  # q flavors CHECK 
ml    = 0.14 # q masses in GeV
el    = [2/3, -1/3, -1/3] # q charge

bk = None  # bk interpolated object

# bk_ is interpolated object
def set_n(bk_):
    global bk
    bk = bk_

# (transverse) wave function for splitting of photon to quark-antiquark dipole
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

def t_integral(z, *args): # args = [qsq2, y]
    m = lambda r_: r_ * psi_t2(z, r_, args[0]) * bk.n(r_, args[1])
    return 2 * np.pi * quad(m, 1e-6, 1e2, epsabs=1e-5, epsrel=0.0)[0]

def l_integral(z, *args): # *args = [qsq2, y]
    m = lambda r_: r_ * psi_l2(z, r_, args[0]) * bk.n(r_, args[1])
    return 2 * np.pi * quad(m, 1e-6, 1e2, epsabs=1e-5, epsrel=0.0)[0]

def t_xsection(x, qsq2, sigma):
    x1 = x * (1 + 4 * ml * ml/qsq2)
    rap   = np.log(x0/x1)
    return 2 * 2.5681 * sigma * quad(t_integral, 0., 0.5, epsabs=1e-4, epsrel=0.0,args=(qsq2, rap))[0]
                
def l_xsection(x, qsq2, sigma):
    x1 = x * (1 + 4 * ml * ml/qsq2)
    rap = np.log(x0/x1)
    return 2 * 2.5681 * sigma * quad(l_integral, 0., 0.5, epsabs=1e-4, epsrel=0.0, args=(qsq2, rap))[0]

def fl(x, qsq2, sigma):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * l_xsection(x, qsq2, sigma)
 
def f2(x, qsq2, sigma):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2, sigma) + l_xsection(x, qsq2, sigma))

def reduced_x(x, qsq2, root_s, sigma):
    y = qsq2/(root_s * root_s * x)  # root_s is center of mass energy
    d = 1 + (1 - y) * (1 - y)

    a = f2(x, qsq2, sigma)
    b = (y * y/d) * fl(x, qsq2, sigma)
    c = a - b
    return [a, b, c]

# saturation scale, small-x values (array type), CM energy sqrt(sNN), interpolated object bk,  observable
def test(q_, x_, cme_, sigma, n_, obs, filename, description=''):
    set_n(n_)

    q      = q_  # GeV ^2
    x      = x_
    sqrt_s = cme_

    f_exist = exists(filename)
    if not f_exist:
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['# bk: ' + n_.get_name()])
            writer.writerow(['q2', 'x', 'cme', 'f2', 'fl', 'redx'])

    with open(filename, 'a') as f:
        writer = csv.writer(f, delimiter='\t')

        if obs == 'f2':
            f2_res = [f2(x[i], q, sigma) for i in range(len(x))]
            for i in range(len(x)):
                writer.writerow([q, x[i], sqrt_s, f2_res[i], '-', '-'])

        elif obs == 'fl':
            fl_res = [fl(x[i], q, sigma) for i in range(len(x))]
            for i in range(len(x)):
                writer.writerow([q, x[i], sqrt_s, '-', fl_res[i], '-'])

        elif obs == 'redx':
            results = [reduced_x(x[i], q, sqrt_s, sigma) for i in range(len(x))]
            f2_res  = [results[i][0] for i in range(len(results))]
            fl_res  = [results[i][1] for i in range(len(results))]
            redx    = [results[i][2] for i in range(len(results))]
            for i in range(len(x)):
                writer.writerow([q, x[i], sqrt_s, f2_res[i], fl_res[i], redx[i]])

'''if __name__ == '__main__':

    with open('params.csv', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        p      = next(reader)
        q      = next(reader)
    
    sig = float(p[0]) * 2 
    bk  = N(p[1], 'dis')
    res = p[2]

    qsq = [float(q[i]) for i in range(len(q))]
    for i in qsq:
        print(i)
        test(i, np.logspace(-5, -2, 20), 319., sig, bk, 'redx', res)'''
