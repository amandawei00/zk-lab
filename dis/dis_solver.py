import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as intg
import scipy.special as spec
import csv
from bk_interpolate_interp2d import N

# Units in GeV
sig = 1.  # constant
alpha = 1./137  # FIND VALUE (EM coupling)
e = 1.60217e-19  # Coulomb
x0 = 0.01  # Bjorken-x, value of x at which evolution starts. (highest experimental value of x included in fit)
q0 = 1.0  # GeV

lamb = 0.241  # lambda_QCD (GeV)
root_sNN = 319 # for data from 1993, s = 4 * Ep * Ee

"""ordering of flavors:
   1. up      -1. antiup
   2. down    -2. antidown
   3. charm   -3. anticharm
   4. strange -4. antistrange
   5. top     -5. antitop
   6. bottom  -6. antibottom  
"""

flavors = [1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6]  # CHECK VALUES (quark flavors, including charm (1.5 GeV))
mf = [0.002, 0.0045, 1.270, 0.101, 172, 5., 0.002, 0.0045, 1.270, 0.101, 172, 5.] # * np.full(1, np.power(3.e8,2)) quark masses in GeV
ef = [2/3, -1/3, 2/3, -1/3, 2/3, -1/3, -2/3, 1/3, -2/3, 1/3, -2/3, 1/3] # * np.full(1, self.e)  # CHECK VALUES quark charges

n = N()  # bk interpolated

# (transverse) wave function for splitting of photon to quark-antiquark dipole
def psi_t2(z, r, *args):
    coeff = (6 * alpha)/(4 * np.pi * np.pi)
    sum = 0

    for i in range(len(flavors)):   # summing over all flavors
        eta2 = eta_squared(z, mf[i], args[0])
        eta = np.sqrt(eta2)
        k02 = np.power(spec.kn(0, eta * r), 2) # modified Bessel function (2nd kind, 0th order)
        k12 = np.power(spec.kn(1, eta * r), 2)  # MacDonald's Function first order

        t1 = (np.power(z, 2) + np.power(1-z, 2)) * eta2 * k12
        t2 = np.power(mf[i], 2) * k02
        sum += np.power(ef[i], 2) * (t1 + t2)
    return coeff * sum

# (longitudinal) wave function for splitting of photon to quark-antiquark dipole
def psi_l2(z, r, *args): # args = [qsq2]

    coeff = (6 * alpha) / (4 * np.power(np.pi, 2))
    sum = 0

    for i in range(len(flavors)):
        eta2 = eta_squared(z, mf[i], args[0])
        eta = np.sqrt(eta2)
        k02 = np.power(spec.kn(0, eta * r), 2)
        sum += 4 * args[0] * np.power(z * (1 - z), 2) * k02 * np.power(ef[i], 2)
    return coeff * sum
        
def eta_squared(z, m_f, qsq2):
    return z * (1 - z) * qsq2 + np.power(m_f, 2)

def t_integral(z, *args):
    m = lambda r_: psi_t2(z, r_, args[0]) * n.master(r_, args[1])
    return intg.quad(m, 0.0, 1/args[0], epsabs=1.e-5)[0]

def l_integral(z, *args): # *args = [qsq2, y]
    m = lambda r_: psi_l2(z, r_, args[0]) * n.master(r_, args[1])
    return intg.quad(m, 0.0, 1/args[0], epsabs=1.e-5)[0]

def t_xsection(x, qsq2):
    r = 1/qsq2
    r0 = (1/q0) * np.power(x/x0, lamb/2)
    y = np.log(x0/x)
    return 2 * np.pi * sig * intg.quad(t_integral, 0., 1., epsabs=1.e-5, args=(qsq2, y))[0]
                
def l_xsection(x, qsq2):
    r = 1/qsq2
    r0 = (1/q0) * np.power(x/x0, lamb/2)  # where does q0 come from?
    y = np.log(x0/x)  # where does y come from? //2 * np.pi comes from angular independence of inner integral
    return 2 * np.pi * sig * intg.quad(l_integral, 0., 1., epsabs=1.e-5, args=(qsq2, y))[0]

def fl(self, x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2) + l_xsection(x, qsq2))
 
def f2(x, qsq2):
    prefac = qsq2/(4 * np.pi * np.pi * alpha)
    return prefac * (t_xsection(x, qsq2) + l_xsection(x, qsq2))


'''
c = 2.568  # unit conversion factor
q = 12  # GeV ^2
x = [0.000261, 0.000383, 0.000562, 0.000825, 0.00133, 0.00237, 0.00421, 0.0075, 0.0133]
s = np.power(296, 2)

f = [f2(x[i], q) for i in range(len(x))]

for i in range(len(x)):
    print("x = " + str(x[i]) + ", f2 = " + str(f[i]))
'''
