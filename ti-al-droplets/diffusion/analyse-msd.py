"""
    Script to analyze MD simulation output in order to compute the shear viscosity coefficient, 
    based on Einstein formula (loosely based on how Gromacs gmx_energy code).
"""

import numpy as np
import matplotlib.pyplot as plt

def read_xvg(file_path):
    t = []
    x = []
    n = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith('@'):
                continue
            n += 1
            values = line.split()
            if not values:
                continue
            t.append(float(values[0]))
            if not x:
                x = [[] for _ in range(len(values) - 1)]
            for i, value in enumerate(values[1:]):
                x[i].append(float(value))
    return t, x, n

num_sets = 3            # Number of MSD components
dt = 0.001              # Time step [ps]

# Reading .xvg output
iframe, msd, n = read_xvg('000/msd-1.dat')
# iframe, msd, n = read_xvg('100/msd-2.dat')

t = [i*dt for i in iframe]

# Technically, the fit should be perfromed only after
# the 'balistic' time, which is usually very short
ti = t[0]   # ps
tf = t[-1]  # ps
idx_ti = int(len(t)*(ti/t[-1]))
idx_tf = int(len(t)*(tf/t[-1]))
tfit = t[idx_ti:idx_tf]

p = np.polyfit(tfit,msd[num_sets][idx_ti:idx_tf],deg=1)
diff_coeff = p[0]/6

print("# ----------------------- #")
print("D =", diff_coeff, "A/ps")
print("# ----------------------- #")

# for m in range(num_sets) :
#     plt.plot(tint, avint[m])
plt.plot(t, msd[num_sets], 'ko')
plt.plot(tfit, np.polyval(p,tfit), 'r-')
plt.xlabel(r'$t$ [ps]')
plt.ylabel(r'$MSD$ [A$^2$]')
plt.show()

# Characteristic length and time
L = 10                  # [A]
tau = L*L/diff_coeff    # [ps]
print("tau =", tau, "ps")