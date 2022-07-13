# 1. set variables in parameters.txt
import sys
import numpy as np
import scipy.integrate as integrate
import csv
from os.path import exists

sys.path.append('python_scripts')
sys.path.append('../bk/')
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
        print(self.flavors)
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
        bkf = self.n.nff(q/z, x2)
        ff_hq = self.ff.xfxQ2(self.f, z, q2)/z

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
        bka = self.n.nfa(q/z, x2)
        ff_hg = self.ff.xfxQ2(self.f, z, q2)/z

        return (1/(z * z)) * pdf_gp * bka * ff_hg

################################################################################################################################
    def rhs(self,pt): # DEFINE TOL
        self.p_t = pt
        xf = self.get_xf()
        m = 0.0

        gluon_intg = integrate.quad(self.integrand1, xf, 1.0, epsabs=1e-5, epsrel=0.0)[0]
        for i in self.flavors:
            if i != 21:
            	self.f = i
            	quark = integrate.quad(self.integrand,xf,1.0, epsabs=1e-5, epsrel=0.0)[0]
            	m += quark + gluon_intg # integral
	
        turkey = m * self.K/(4 * np.pi * np.pi)
        print("rhs exit, pt = " + str(pt) + ", rhs = " + str(turkey))

        return turkey
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if __name__=='__main__':

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

    pro  = params[0]
    h    = params[1]
    y    = float(params[2])
    qsq2 = float(params[3])
    snn  = float(params[4])
    k    = float(params[5])
    a    = float(params[6])
    b    = float(params[7])
    n    = int(params[8])
    init = params[9]
    res  = params[10]
    orde = params[11]

    from_file = '../bk/results/' + 'RK' + orde + '/' + init
    to_file   = 'results/' + pro + '/' + 'RK' + orde + '/' + res

    s = Master(h, y, qsq2, snn, k, from_file)  # creating instance of class
    dp_t = (b - a)/n

    p_t = np.arange(a, b, dp_t)  # tranverse momenta values
    cs = np.zeros(len(p_t))

    if not exists(to_file):
        with open(to_file, 'w') as tfile: # write temporary output file
            writer = csv.writer(tfile, delimiter='\t')
            writer.writerow(['# ', 'par', 'q2', 'pt_low', 'pt_high', 'order'])
            writer.writerow(['# ', h, qsq2, a, b, orde])
            writer.writerow(['cme', 'y', 'pt', 'dN'])
            for i in range(len(p_t)):
                cs[i] = s.rhs(p_t[i])
                writer.writerow([snn, y, p_t[i], cs[i]])
    else:
        with open(to_file, 'a') as tfile:
            writer = csv.writer(tfile, delimiter='\t')
            for i in range(len(p_t)):
                cs[i] = s.rhs(p_t[i])
                writer.writerow([snn, y, p_t[i], cs[i]])

# end of program
