from scipy import LowLevelCallable as llc
import scipy.integrate as intg
# import test
import t
import subprocess
import numpy as np

# subprocess.run(['python3', 'setup.py', 'build_ext', '--inplac'])

def sample(x, y):
    return x * x + 2 * y

single_intg = llc.from_cython(t, 'py_f3')
print(intg.quad(single_intg, 0, 4))

print("python integral")
print(intg.nquad(sample, [[0, 4], [0, 2]]))

double_intg = llc.from_cython(t, 'py_f2', signature='double (int, double *)')
print("cython integral")

# x * x + 2 * y, xx=[x,y]. assuming x is first
print(intg.dblquad(double_intg, 0, 6, 0, 2, epsabs=0.0, epsrel=0.5))
print(intg.dblquad(double_intg, 0, 2, 0, 6, epsabs=0.0, epsrel=0.5))
#print(intg.nquad(double_intg, [[0, 4],[0,2]]))

# print(t.py_f(np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.intc), 7))
# t.set_abc(3,2, 5)
# t.print_abc()
# t.set_abc(2, 2, 2)
# t.print_abc()
