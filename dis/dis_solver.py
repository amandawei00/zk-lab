import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.special as spec
import csv

from os.path import exists

sys.path.append('../bk/')
from bk_interpolate import N

# Units in GeV
alpha = 1./137   # FIND VALUE (EM coupling)
x0    = 0.01
lamb  = 0.241  # lambda_QCD (GeV)

"""ordering of flavors:
   1. up      -1. antiup
   2. down    -2. antidown
   3. strange -3. strange
   4. charm   -4. anticharm
   5. top     -5. antitop
   6. bottom  -6. antibottom  
"""

light = [1, 2, 3, -1, -2, -3]  # q flavors CHECK 
ml    = [0.002, 0.0045, 1.270, 0.002, 0.0045, 1.270] # q masses in GeV
el    = [2/3, -1/3, -1/3, -2/3, 1/3, 1/3] # q charge

heavy = [4, 5, 6, -4, -5, -6]
mh    = [1.270, 172., 5., 1.270, 172., 5.]
eh    = [2/3, 2/3, -1/3, -2/3, -2/3, 1/3]

bk = None  # bk interpolated object

# n_ is interpolated object
def set_n(bk_):
    global bk
    bk = bk_

# (transverse) wave function for splitting of photon to quark-antiquark dipole
def psi_t2(z, r, *args):
    coeff = (6 * alpha)/(4 * np.pi * np.pi)
    sum   = 0

    for i in range(len(light)):   # summing over all flavors
        eta2 = eta_squared(z, ml[i], args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(spec.kn(0, eta * r), 2) # modified Bessel function (2nd kind, 0th order)
        k12  = np.power(spec.kn(1, eta * r), 2)  # MacDonald's Function first order

        t1   = (z * z + (1 - z) * (1 - z)) * eta2 * k12
        t2   = ml[i] * ml[i] * k02
        sum  += el[i] * el[i] * (t1 + t2)
    return coeff * sum

# (longitudinal) wave function for splitting of photon to quark-antiquark dipole
def psi_l2(z, r, *args): # args = [qsq2]

    coeff = (6 * alpha) / (4 * np.pi * np.pi)
    sum   = 0

    for i in range(len(light)):
        eta2 = eta_squared(z, ml[i], args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(spec.kn(0, eta * r), 2)
        sum  += 4 * args[0] * np.power(z * (1 - z), 2) * k02 * el[i] * el[i]
    return coeff * sum
        
def eta_squared(z, m_f, qsq2):
    return z * (1 - z) * qsq2 + m_f * m_f

def t_integral(z, *args):
    m = lambda r_: r_ * psi_t2(z, r_, args[0]) * bk.n(r_, args[1])
    return quad(m, 3.e-6, 60., epsabs=1.e-3, epsrel=0.0)[0]

# orignal integration bound: [3.e-6, 1/args[0]]
def l_integral(z, *args): # *args = [qsq2, y]
    m = lambda r_: r_ * psi_l2(z, r_, args[0]) * bk.n(r_, args[1])
    return quad(m, 3.e-6, 60., epsabs=1.e-3, epsrel=0.0)[0]

def t_xsection(x, qsq2, sigma):
    y   = np.log(x0/x)
    return 2 * np.pi * sigma * quad(t_integral, 0., 1., epsabs=1.e-3, epsrel=0.0,  args=(qsq2, y))[0]
                
def l_xsection(x, qsq2, sigma):
    y  = np.log(x0/x)  # 2 * np.pi comes from angular independence of inner integral
    return 2 * np.pi * sigma * quad(l_integral, 0., 1., epsabs=1.e-3, epsrel=0.0,  args=(qsq2, y))[0]

def fl(x, qsq2, sigma):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * l_xsection(x, qsq2, sigma)
 
def f2(x, qsq2, sigma):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2, sigma) + l_xsection(x, qsq2, sigma))

# at low qsq (qsq << MZ^2, and Z exchange is negligible)
def reduced_x(x, qsq2, root_s, sigma):
    y = qsq2/(root_s * root_s * x)  # root_s is center of mass energy
    d = 1 + (1 - y) * (1 - y)

    a = f2(x, qsq2, sigma)
    b = (y * y/d) * fl(x, qsq2, sigma)
    c = a - b
    return [a, b, c]

# saturation scale, small-x values (array type), CM energy sqrt(sNN), interpolated object bk,  observable
# c = 2.58 unit conversion factor?
def test(q_, x_, cme_, sigma, n_, obs, filename, description=''):
    set_n(n_)

    q      = q_  # GeV ^2
    x      = x_
    sqrt_s = cme_

    f_exist = exists(filename)
    if not f_exist:
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([description])
            writer.writerow(['q2', 'x', 'f2', 'fl', 'redx'])

    with open(filename, 'a') as f:
        writer = csv.writer(f, delimiter='\t')

        if obs == 'f2':
            f2_res = [f2(x[i], q) for i in range(len(x))]
            for i in range(len(x)):
                writer.writerow([q, x[i], f2_res[i], '-', '-'])

        elif obs == 'fl':
            fl_res = [fl(x[i], q) for i in range(len(x))]
            for i in range(len(x)):
                writer.writerow([q, x[i], '-', fl_res[i], '-'])

        elif obs == 'redx':
            results = [reduced_x(x[i], q, sqrt_s, sigma) for i in range(len(x))]
            f2_res  = [results[i][0] for i in range(len(results))]
            fl_res  = [results[i][1] for i in range(len(results))]
            redx    = [results[i][2] for i in range(len(results))]
            for i in range(len(x)):
                writer.writerow([q, x[i], f2_res[i], fl_res[i], redx[i]])
                print("x = " + str(x[i]) + ", redx = " + str(redx[i]))

# bk = N('../bk/results/bk_MV1.csv') 
# test(45.0, np.logspace(-6, -2, 25), 318.0, 28, bk, 'redx', 'test.csv')
