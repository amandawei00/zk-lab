from numpy import exp, log

# qsq0 [GeV^2] = 0.6 [GeV^2]: initial saturation scale
# x0           = 0.01       : onset of small-x
# gamma        = 0.3        : 
def s(x, r_, x0=0.01, qsq0=1.0 gamma=0.3):
    qsq = qsq0 * np.power(x0/x, gamma)
    return exp(-(1/4) * qsq * r_ * r_ * log(1/(gamma * r_) + exp(1))
