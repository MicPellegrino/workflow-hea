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

folder = '000'
n_rep = 9

# Reading and averaging .xvg output
iframe, vac_dummy, n = read_xvg(folder+'/vacf-1-0.dat')
iframe, ivac_dummy, n = read_xvg(folder+'/ivacf-1-0.dat')
vac1_m = np.zeros_like(vac_dummy)
vac2_m = np.zeros_like(vac_dummy)
ivac1_m = np.zeros_like(ivac_dummy)
ivac2_m = np.zeros_like(ivac_dummy)
for nr in range(n_rep) :
    frame, vac1, n = read_xvg(folder+'/vacf-1-'+str(nr)+'.dat')
    iframe, ivac1, n = read_xvg(folder+'/ivacf-1-'+str(nr)+'.dat')
    frame, vac2, n = read_xvg(folder+'/vacf-1-'+str(nr)+'.dat')
    iframe, ivac2, n = read_xvg(folder+'/ivacf-1-'+str(nr)+'.dat')
    vac1_m += vac1
    ivac1_m += ivac1
    vac2_m += vac2
    ivac2_m += ivac2

vac1_m /= n_rep
vac2_m /= n_rep

vac1_m = vac1_m / vac1_m[:,0][:,None]
vac2_m = vac2_m / vac2_m[:,0][:,None]

ivac1_m /= n_rep
ivac2_m /= n_rep

t = [i*dt for i in frame]
it = [i*dt for i in iframe]

fig, (ax1,ax2) = plt.subplots(1,2)

ax1.plot(t, vac1_m[num_sets])
ax1.set_xlabel(r'$t$ [ps]')
ax1.set_ylabel(r'$norm. VACF$ []')
ax1.set_xlim([t[0],t[-1]])

ax2.plot(it, ivac1_m[0])
ax2.set_xlabel(r'$t$ [ps]')
ax2.set_ylabel(r'$I(VACF)$ [A$^2$/ps]')
ax2.set_xlim([it[0],it[-1]])

plt.show()