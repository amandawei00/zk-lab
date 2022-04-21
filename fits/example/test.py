import sys,os
import numpy as np
import pylab as py
import pandas as pd
from scipy.integrate import quad,fixed_quad
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib
from matplotlib.pyplot import gca
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
from scipy.interpolate import interp1d
import random
sizeOfFont = 20
rc('text',usetex=True)
fontProperties = {'weight' : 'normal', 'size' : sizeOfFont}
rc('text',usetex=True)
rc('font',**fontProperties)
from scipy.interpolate import interp1d
from iminuit import Minuit
from iminuit.cost import LeastSquares
import warnings
warnings.filterwarnings("ignore")

dfpythia = pd.read_csv("example.dat")
midpoints = np.array(dfpythia.jT)
errors = np.array(dfpythia.sigma)
values = np.array(dfpythia.exp)

# construct theory
def toy_thy(pT,Nn,pT2):
    return Nn*pT*np.exp(-pT*pT/pT2)

# perform fit
def chi_squared(Nn,pT2):
    res = 0.
    for i in range(len(midpoints)):
        jT = midpoints[i]
        if jT <= 2.:
            theory = toy_thy(jT,Nn,pT2)
            pythia = values[i]
            errval = errors[i]
            res += (theory-pythia)**2./errval**2.0
    return res
chi_squared.errordef = Minuit.LEAST_SQUARES
m = Minuit(chi_squared,Nn = 0.,pT2= 0.)
m.simplex()
print(m.values)  # prints fit results
print(m.errors)  # prints errors
# print(m.fval/(len(midpoints) - len(m.values)))  # prints reduced chi2 (goodness of fit)
print(repr(m.fmin))
