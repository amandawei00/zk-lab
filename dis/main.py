import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as intg
import scipy.special as spec
import csv
from bk_interpolate_interp2d import N

# Units chosen to be in GeV


class Solve:

    def __init__(self, tl):
        self.sig = 1.
        self.alpha = 1./137  # FIND VALUE (EM coupling)
        self.e = 1.60217e-19  # Coulomb
        self.x0 = 0.01  # value of x at which evolution starts. (highest experimental value of x included in fit)

        """ordering of flavors:
            1. up      -1. antiup
            2. down    -2. antidown
            3. charm   -3. anticharm
            4. strange -4. antistrange
            5. top     -5. antitop
            6. bottom  -6. antibottom  
        """

        self.f = [1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6]  # CHECK VALUES (quark flavors, including charm with mass of 1.5 GeV)
        self.mf = [0.002, 0.0045, 1.270, 0.101, 172, 5., 0.002, 0.0045, 1.270, 0.101, 172, 5.] # * np.full(1, np.power(3.e8,2)) quark masses in GeV
        self.ef = [2/3, -1/3, 2/3, -1/3, 2/3, -1/3, -2/3, 1/3, -2/3, 1/3, -2/3, 1/3] * np.full(1, self.e)  # CHECK VALUES quark charges
        self.tl = tl

        self.n = N()  # bk interpolated

        # these variables are initialized here. to be assigned in self.rhs()
        self.x = 0.0
        self.qsq2 = 0.0
        self.y = 0.0  # rapidity

        self.r0 = 0.0

    # (transverse) wave function for splitting of photon to quark-antiquark dipole
    def psi_t(self, z, r):
        coeff = (6 * self.alpha)/(4 * np.power(np.pi, 2))
        sum = 0

        for i in range(len(self.f)):   # summing over all flavors
            eta2 = self.eta_squared(z, self.mf[i])
            eta = np.power(eta2, 0.5)
            k0 = spec.kn(0, eta * r)  # modified Bessel function of the second kind, 0th order (MacDonald's Function)
            k1 = spec.kn(1, eta * r)  # MacDonald's Function first order

            t1 = (np.power(z, 2) + np.power(1-z, 2)) * eta2 * np.power(k1, 2)
            t2 = np.power(self.mf[i], 2) * np.power(k0, 2)

            sum += np.power(self.ef[i], 2) * (t1 + t2)

            # print("eta2 = " + str(eta2) + ", k0 = " + str(k0) + ", k1 = " + str(k1) + ", t1 = " + str(t1) + ", t2 = " + str(t2))

        # print("psi_t for z = " + str(z) + ", r = " + str(r) + " is: " + str(coeff * sum))
        return coeff * sum

    # (longitudinal) wave function for splitting of photon to quark-antiquark dipole
    def psi_l(self, z, r):

        coeff = (6 * self.alpha) / (4 * np.power(np.pi, 2))
        sum = 0

        for i in range(len(self.f)):
            eta2 = self.eta_squared(z, self.mf[i])
            eta = np.power(eta2, 0.5)
            k0 = spec.kn(0, eta * r)
            sum += np.power(2 * self.ef[i] * self.qsq2 * z * (1 - z) * k0, 2)
        
        return coeff * sum
        
    def eta_squared(self, z, m_f):
        # print("z = " + str(z) + ", qsq2 = " + str(qsq2) + ", m_f = " + str(m_f))
        return z * (1 - z) * self.qsq2 + np.power(m_f, 2)

    def inner_integral(self, z):
        if self.tl == "T":
            m = lambda r_: np.power(np.abs(self.psi_t(z, r_)), 2) * self.n.master(r_, self.y)
        elif self.tl == "L":
            m = lambda r_: np.power(np.abs(self.psi_l(z, r_)), 2) * self.n.master(r_, self.y)

        """ t = np.linspace(0, 1000)
        yt = m(t)
        plt.xlim(0, 1000)
        plt.ylim(0, 1)
        plt.plot(t, yt)
        plt.show() """

        u = intg.quad(m, 0.0, self.r0, epsabs=1.e-5)[0]
        print("inner integral for z = " + str(z) + " u is: " + str(u))
        return u

    def rhs(self, x, qsq2):
        self.x = x
        self.qsq2 = qsq2
        self.r0 = np.sqrt(1/qsq2)

        self.y = np.log(self.x0/self.x)
        integral = intg.quad(self.inner_integral, 0, 1, epsabs=1.e-5)[0]
        return 2 * np.pi * self.sig*integral  # 2*np.pi comes from theta component in the inner double integral. since there is no dependence on theta, multiply it out


if __name__ == "__main__":
    # p = []  # p for parameters

    """ order of parameters in parameters.txt:
        1. transverse/longitudinal momentum """
    """# reading in values from parameters.txt
    with open("parameters.txt", "r") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for r in reader:
            raw = r[0].split()
            p.append(raw[len(raw)-1])"""

    # create new instance of solve class
    s = Solve("L")
    # p = Solve("L")
    # x = np.logspace(10.e-3, 10.e-3, num=4)
    x = [10.e-3]
    qsq2 = 10.0  # GeV^2
    for i in range(len(x)):
        a = s.rhs(x[i], qsq2)
        print("for x = " + str(x[i]) + ", sigma = " + str(a))