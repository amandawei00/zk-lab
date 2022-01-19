import matplotlib.pyplot as plt
from scipy.integrate import simps
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

''' template
qsq2 = 
x_th = []
y_th = []

x_exp = x_th[:]
y_exp = []

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/2.7 for i in range(len(y_th))]

exp_error_stat = []
exp_error_sys = []

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]

'''

def normalize(x_th, y_th, x_exp, y_exp):
    f_th = InterpolatedUnivariateSpline(x_th, y_th, k=1)
    f_exp = InterpolatedUnivariateSpline(x_exp, y_exp, k=2)

    a_th = f_th.integral(x_th[0], x_th[len(x_th)-1])
    a_exp = f_exp.integral(x_exp[0], x_exp[len(x_exp)-1])

    return float(a_exp/a_th)


#qsq2 = 2.5 GeV^2
"""
x_th = [0.0001, 0.000177827941, 0.000316227766, 0.0005623413252, 0.001, 0.00177827941, 0.00316227766, 0.005623413252, 0.01]
y_th = [1.886886362, 1.704623514, 1.537912721, 1.385617549, 1.246655295, 1.12002439, 1.00475126, 0.9000246767, 0.8048919988]

y_th = [y_th[i] * (1./1.88) for i in range(len(y_th))]
"""


# qsq2 = 12 GeV^2
'''
qsq2 = 12
x_th = [0.000261, 0.000383, 0.000562, 0.000825, 0.00133, 0.00237, 0.00421, 0.00750, 0.01330]
y_th = [3.872595062, 3.516630264, 3.186674368, 2.880772742, 2.532035557, 2.153612931, 1.820225504, 1.524204608, 1.390125278]

x_exp = x_th[:]
y_exp = [1.35, 1.26, 1.19, 1.08, 0.96, 0.85, 0.74, 0.70, 0.58]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/2.7 for i in range(len(y_th))]

exp_error_stat = [0.06, 0.05, 0.05, 0.04, 0.05, 0.04, 0.04, 0.04, 0.04]
exp_error_sys = [0.13, 0.12, 0.12, 0.11, 0.13, 0.1, 0.1, 0.11, 0.12]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 20 GeV^2
'''
qsq2 = 20
x_th = [0.000562, 0.000825, 0.00133, 0.00237, 0.00421, 0.00750, 0.01330, 0.0237]
y_th = [3.897583272, 3.492729625,3.034960216,2.543476252,2.115846024,1.741236143,1.573421918,1.573421918]

x_exp = x_th[:]
y_exp = [1.52, 1.17, 1.03, 1.03, 0.83, 0.74, 0.64, 0.51]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/2.7 for i in range(len(y_th))]

exp_error_stat = [0.08, 0.07, 0.05, 0.05, 0.04, 0.04, 0.04, 0.05]
exp_error_sys = [0.12, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.08]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 50 GeV^2
'''
qsq2 = 50.
x_th = [0.00133, 0.00237, 0.00421, 0.00750, 0.01330, 0.0237, 0.0422]
y_th = [ 4.141828944,
3.400898628,
2.765487118,
2.221071047,
1.981018843,
1.981018843,
1.981018843]

x_exp = x_th[:]
y_exp = [1.46, 1.08, 1.00, 0.65, 0.66, 0.52, 0.40]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/3.0 for i in range(len(y_th))]

exp_error_stat = [0.10, 0.08, 0.07, 0.06, 0.06, 0.05, 0.05]
exp_error_sys = [0.12, 0.09, 0.09, 0.08, 0.08, 0.07, 0.08]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 65 GeV^2
'''
qsq2 = 65.
x_th = [0.00237, 0.00421, 0.0075, 0.0133]
y_th = [3.703571381,
2.998123702,
2.393803854,
2.129127745]

x_exp = x_th[:]
y_exp = [1.4, 1.09, 0.95, 0.69]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/2.7 for i in range(len(y_th))]

exp_error_stat = [0.11, 0.09, 0.08, 0.08]
exp_error_sys = [0.11, 0.09, 0.11, 0.07]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 80 GeV^2
'''
qsq2 = 80.
x_th = [0.00237, 0.00421, 0.0075, 0.0133]
y_th = [3.961462275,
3.19465731,
2.54138463,
2.255435131]

x_exp = x_th[:]
y_exp = [1.19, 1.09, 0.7, 0.71]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/3.2 for i in range(len(y_th))]

exp_error_stat = [0.13, 0.14, 0.08, 0.07]
exp_error_sys = [0.13, 0.12, 0.06, 0.08]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 120 GeV^2
'''
qsq2 = 120
x_th = [0.00237, 0.00421, 0.00750, 0.01330]
y_th = [4.496041345,
3.600790834,
2.841658065,
2.513390747]

x_exp = x_th[:]
y_exp = [1.60, 0.99, 0.83, 0.73]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/3.2 for i in range(len(y_th))]

exp_error_stat = [0.21, 0.14, 0.12, 0.11]
exp_error_sys = [0.15, 0.14, 0.13, 0.11]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# qsq2 = 200 GeV^2
'''
qsq2 = 200.
x_th = [0.00421, 0.0075, 0.0133]
y_th = [4.075703402,
3.196068555,
2.813441222]

x_exp = x_th[:]
y_exp = [1.41, 0.91, 0.72]

a = normalize(x_th, y_th, x_exp, y_exp)
y_th = [y_th[i]/3.2 for i in range(len(y_th))]

exp_error_stat = [0.14, 0.1, 0.09]
exp_error_sys = [0.11, 0.11, 0.08]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]
'''

# Q^2 = 300 GeV^2 (319 GeV)
qsq2=300
x_th = [0.008, 0.013, 0.02, 0.032, 0.05, 0.08, 0.13, 0.18, 0.25, 0.4]
y_th = [3.652134081, 3.393053954, 3.260312703, 3.166112403, 3.094053085, 3.050276371, 3.028754695, 3.021628475, 3.017566601, 3.014637041]

x_exp = x_th
y_exp = [0.9921, 0.8204, 0.6985, 0.5863, 0.4965, 0.4188, 0.3551, 0.2953, 0.2804, 0.1477]

y_th = [normalize(x_th, y_th, x_exp, y_exp) * y_th[i] for i in range(len(y_th))]


# qsq2 = 400 GeV^2 (319 GeV)

qsq2 = 400.
x_th = [0.0075, 0.0133]
y_th = [3.576329732, 3.131275517]

x_exp = x_th[:]
y_exp = [1.16, 0.81]

y_th = [y_th[i]/3.3 for i in range(len(y_th))]

exp_error_stat = [0.17, 0.11]
exp_error_sys = [0.08, 0.08]

total_error = [np.sqrt(np.power(exp_error_sys[i],2) + np.power(exp_error_stat[i],2)) for i in range(len(exp_error_stat))]


# qsq2 = 8000 GeV^2 (319 GeV)
"""
x_th = [0.18, 0.25, 0.4, 0.65]
y_th = [4.904070283, 4.61380915, 4.347081647, 4.189467577]

x_exp = x_th[:]
y_exp = [0.2781, 0.2022, 0.1006, 0.0103]

a = normalize(x_th, y_th, x_exp, y_exp)
#y_th = [y_th[i] * a for i in range(len(y_th))]
y_th = [y_th[i]-4 for i in range(len(y_th))]

############################################################################################
"""

############################################################################################

# limits of x-axis
x1 = 0.006
x2 = 0.015

# limits of y-axis
y1 = 0.
y2 = 1.6

# scale of axes
x_scale = 'log'
y_scale = 'linear'

# design parameters
marker_th = 'v'
marker_exp = 'D'
linewidth_th = 1.5

linestyles = ['', '-', '--', '-.', ':']
linestyle_th = linestyles[2]
linestyle_exp = linestyles[0]

x_label = "x"
y_label = "F2(x, Q^2)"

############################################################################################

plt.xlim(x1, x2)
plt.ylim(y1, y2)

plt.xscale(x_scale)
plt.yscale(y_scale)

plt.xlabel("x", fontsize=16)
plt.ylabel("F2", fontsize=16)

plt.title("Q^2 = " + str(qsq2) + " GeV^2", fontsize=16)

plt.plot(x_th, y_th, marker=marker_th, linestyle=linestyle_th, color='magenta', label="calculated")
plt.plot(x_exp, y_exp, marker=marker_exp, linestyle=linestyle_exp, color='blue', label="experimental")

# plt.errorbar(x_exp, y_exp, yerr=total_error, linestyle=linestyle_exp, color='blue', capsize=2)
plt.grid(alpha=.5, linestyle=':')

plt.legend()
plt.show()

