import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as intg
import scipy.special as spec
import csv
from bk_interpolate import N

# Units in GeV
alpha = 1./137   # FIND VALUE (EM coupling)
x0    = 0.01
sig   = 1.
lamb  = 0.241  # lambda_QCD (GeV)

"""ordering of flavors:
   1. up      -1. antiup
   2. down    -2. antidown
   3. charm   -3. anticharm
   4. strange -4. antistrange
   5. top     -5. antitop
   6. bottom  -6. antibottom  
"""

flavors = [1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6]  # q flavors CHECK 
mf      = [0.002, 0.0045, 1.270, 0.101, 172, 5., 0.002, 0.0045, 1.270, 0.101, 172, 5.] # q masses in GeV
ef      = [2/3, -1/3, 2/3, -1/3, 2/3, -1/3, -2/3, 1/3, -2/3, 1/3, -2/3, 1/3] # q charge

n = N('../bk/results/fit1.csv')  # bk interpolated

def set_param(q, s_):
    global qsq0, sig
    qsq0  = q
    sig = s_

# (transverse) wave function for splitting of photon to quark-antiquark dipole
def psi_t2(z, r, *args):
    coeff = (6 * alpha)/(4 * np.pi * np.pi)
    sum   = 0

    for i in range(len(flavors)):   # summing over all flavors
        eta2 = eta_squared(z, mf[i], args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(spec.kn(0, eta * r), 2) # modified Bessel function (2nd kind, 0th order)
        k12  = np.power(spec.kn(1, eta * r), 2)  # MacDonald's Function first order

        t1   = (z * z + (1 - z) * (1 - z)) * eta2 * k12
        t2   = mf[i] * mf[i] * k02
        sum  += ef[i] * ef[i] * (t1 + t2)
    return coeff * sum

# (longitudinal) wave function for splitting of photon to quark-antiquark dipole
def psi_l2(z, r, *args): # args = [qsq2]

    coeff = (6 * alpha) / (4 * np.pi * np.pi)
    sum   = 0

    for i in range(len(flavors)):
        eta2 = eta_squared(z, mf[i], args[0])
        eta  = np.sqrt(eta2)
        k02  = np.power(spec.kn(0, eta * r), 2)
        sum  += 4 * args[0] * np.power(z * (1 - z), 2) * k02 * ef[i] * ef[i]
    return coeff * sum
        
def eta_squared(z, m_f, qsq2):
    return z * (1 - z) * qsq2 + m_f * m_f

def t_integral(z, *args):
    m = lambda r_: r_ * psi_t2(z, r_, args[0]) * n.master(r_, args[1])
    return intg.quad(m, 3.e-6, 60., epsabs=1.e-3)[0]

# orignal integration bound: [3.e-6, 1/args[0]]
def l_integral(z, *args): # *args = [qsq2, y]
    m = lambda r_: r_ * psi_l2(z, r_, args[0]) * n.master(r_, args[1])
    return intg.quad(m, 3.e-6, 60., epsabs=1.e-3)[0]

def t_xsection(x, qsq2):
    y = np.log(x0/x)
    return 2 * np.pi * sig * intg.quad(t_integral, 0., 1., epsabs=1.e-3, args=(qsq2, y))[0]
                
def l_xsection(x, qsq2):
    y  = np.log(x0/x)  # where does y come from? //2 * np.pi comes from angular independence of inner integral
    return 2 * np.pi * sig * intg.quad(l_integral, 0., 1., epsabs=1.e-3, args=(qsq2, y))[0]

def fl(x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * l_xsection(x, qsq2)
 
def f2(x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2) + l_xsection(x, qsq2))

# at low qsq (qsq << MZ^2, and Z exchange is negligible)
def reduced_x(x, qsq2, root_s):
    y = qsq2/(root_s * root_s * x)  # root_s is center of mass energy
    d = 1 + (1 - y) * (1 - y)

    a = f2(x, qsq2)
    b = (y * y/d) * fl(x, qsq2)
    c = a - b
    return [a, b, c]


# c = 2.568  # unit conversion factor
q = 120  # GeV ^2
x = np.logspace(-5, -2, 25)
sqrt_s = 296

f2_res = [f2(x[i], q) for i in range(len(x))]
# f2_res  = [results[i][0] for i in range(len(results))]
# fl_res  = [results[i][1] for i in range(len(results))]
# redx    = [results[i][2] for i in range(len(results))]

for i in range(len(x)):
    print("x = " + str(x[i]) + ", sig = " + str(f2_res[i]))

with open('fit1_results.csv', 'a') as f:
    writer= csv.writer(f, delimiter='\t')
    # writer.writerow(['# comparison of f2 results from bk solution with parameters from fit1.'])
    for i in range(len(x)):
        writer.writerow([sqrt_s, q, x[i], f2_res[i]])
