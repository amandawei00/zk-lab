# 1. set variables in parameters.txt
# 2. verify pdf, ff, N are correct
# 3. check destination filename
# 4. TODO: MAKE 2 and 3 PART OF 1

import sys
import numpy as np
import scipy.integrate as integrate
import csv
# import matplotlib.pyplot as plt
# import subprocess

sys.path.append("python_scripts")
import lhapdf as pdf
from bk_interpolate import N

# - IH FOR HADRON TYPE, IC FOR HADRON CHARGE SHOULD BE MODIFIABLE AND INITIATED UPON CONSTRUCTION

class Master():
    def __init__(self, h, y, qsq, s_NN_root, K, initials):
        self.p = pdf.mkPDF("CT10",0)
        # self.p = pdf.mkPDF("CT10/0")

        self.n = N(initials)
        self.ff = pdf.mkPDF("DSS07PI",0)

        self.qsq2 = qsq # saturation scale
        self.sNN_root = s_NN_root # collision energy per nucleon [GeV]
        self.K = K
        self.hadron = h

        # self.flavors = [1,2,3,4,5,21] # flavors: d(1), u(2), s(3), c(4), b(5), g(21)
        self.flavors = self.p.flavors()
        self.f = 0.0

        self.p_t = 0.0 # plotting points w respect to p_t
        self.yh = y

#############################################################################################################
    def get_xf(self):
        xf = (self.p_t/self.sNN_root)*np.exp(self.yh)
        return xf
#############################################################################################################
    def integrand(self,z):
        xf = self.get_xf()

        x1 = xf/z
        x2 = x1*np.exp(-2*self.yh)
        q = self.p_t
        q2 = q*q

        pdf_qp = self.p.xfxQ2(self.f, x1, q2) # returns x1*f(x1,pt^2) where f is pdf
        bkf = self.n.udg_f(x2, q/z)
        ff_hq = self.ff.xfxQ2(self.f, z, q2)/z

        # return (1/z) * pdf_qp * bkf * ff_hq
        return (1/(z * z)) * pdf_qp * bkf * ff_hq
###################################################################################################################################
    def integrand1(self, z):
        xf = self.get_xf()

        x1 = xf/z
        x2 = x1*np.exp(-2*self.yh)
        q = self.p_t
        q2 = q*q

        self.f = 21
        pdf_gp = self.p.xfxQ2(self.f, x1, q2)
        bka = self.n.udg_a(x2, q/z)
        ff_hg = self.ff.xfxQ2(self.f, z, q2)/z

        # return (1/z) * pdf_gp * bka * ff_hg
        return (1/(z * z)) * pdf_gp * bka * ff_hg

################################################################################################################################
    def rhs(self,pt): # DEFINE TOL
        self.p_t = pt
        xf = self.get_xf()
        m = 0.0

        gluon_intg = integrate.quad(self.integrand1, xf, 1.0, epsabs=0.0, epsrel=0.05)[0]
        for i in self.flavors:
            if i != 21:
            	self.f = i
            	quark = integrate.quad(self.integrand,xf,1.0, epsabs=0.0, epsrel=0.05)[0]
            	m += quark + gluon_intg # integral
	
        turkey = m * self.K/(4 * np.pi * np.pi)
        print("rhs exit, pt = " + str(pt) + ", rhs = " + str(turkey))

        return turkey
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if __name__=="__main__":

    params = []
    """ order of parameters in parameters.txt:
    1.  hadron type
    2.  y
    3.  qsq2
    4.  s_NN
    5.  K
    6.  p_t (lower bound)
    7.  p_t (upper bound)
    8.  number of points evaluated 
    9.  filename initial conditions
    10. file to write results to"""

    with open('params.csv', 'r') as infile:  # opening parameters file to read in input
        reader  = csv.reader(infile, delimiter='\t')	
        header  = next(reader)
        params  = next(reader)

    h    = params[0]
    y    = float(params[1])
    qsq2 = float(params[2])
    snn  = float(params[3])
    k    = float(params[4])
    a    = float(params[5])
    b    = float(params[6])
    n    = int(params[7])
    init = params[8]
    res  = params[9]

    s = Master(h, y, qsq2, snn, k, init)  # creating instance of class
    dp_t = (b - a)/n

    p_t = np.arange(a, b, dp_t)  # tranverse momenta values
    cs = np.zeros(len(p_t))

    # if res exists:
    # else: 
    with open(res, 'a') as tfile: # write temporary output file
        writer = csv.writer(tfile, delimiter='\t')
        for i in range(len(p_t)):
            cs[i] = s.rhs(p_t[i])
            writer.writerow([float(param[3]), float(param[1]), p_t[i], cs[i]])

# end of program
