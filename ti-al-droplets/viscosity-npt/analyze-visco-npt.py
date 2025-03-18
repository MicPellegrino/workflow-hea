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

diag_fac = 0.75         # Factor for diagonal components
num_sets = 6            # Number of stress components
dt = 0.001              # Time step [ps]
T = 2000                # Temperature [K]
kB = (1e3)*138.0649     # kB [bar*nm^3/K]

# Reading .xvg output
t, eint, nint = read_xvg('100/eint.dat')
t, vol, nint = read_xvg('100/volume.dat')

vol = np.array(vol[0])
V = (1e-3)*np.mean(vol)     # Volume [nm^3]
std_V = (1e-3)*np.std(vol)

scale = V/(2.0*T*kB)

# Number of averaging blocks and restrarts for computing averages
num_blocks = 10
m = 10
num_restarts = nint//(m*num_blocks)
nint_block = nint//num_blocks+1
step_size = max(nint_block//num_restarts,1)

tint = np.zeros(num_restarts+1)
avint = np.zeros((len(eint),num_restarts+1))

print("num_blocks =", num_blocks)
print("num_restarts =", num_restarts)
print("nint_block =", nint_block)
print("step_size =", step_size)

for i in range(0, nint_block, step_size) :
    ii = i//step_size
    print("Averaging step", ii, "/", num_restarts)
    # Off-diagonal components
    for m in range(3) :
        for j in range(0,nint-i) :
            di = (eint[m][j + i]-eint[m][j])**2
            avint[m][ii] += di
            avint[num_sets][ii] += di/num_sets
    # Diagonal components
    for m in range(3,num_sets) :
        for j in range(0,nint-i) :
            di = diag_fac*(eint[m][j + i]-eint[m][j])**2
            avint[m][ii] += di
            avint[num_sets][ii] += di/num_sets
    tint[ii] = dt*t[i]

std_avint = np.std(avint[0:num_sets], axis=0)
std_avint_p = std_avint/avint[num_sets]

# Technically, the fit should be perfromed only between
# the decorrelation time (ti) and the time of max prescribed
# standard error (tf)
ti = 500        # ps
tf = tint[-1]   # ps
idx_ti = int(len(tint)*(ti/tint[-1]))
idx_tf = int(len(tint)*(tf/tint[-1]))
tfit = tint[idx_ti:idx_tf]

p = np.polyfit(tfit,avint[num_sets][idx_ti:idx_tf],deg=1)
p_low = np.polyfit(tfit,(1-std_avint_p[idx_ti:idx_tf])*avint[num_sets][idx_ti:idx_tf],deg=1)
p_upp = np.polyfit(tfit,(1+std_avint_p[idx_ti:idx_tf])*avint[num_sets][idx_ti:idx_tf],deg=1)
viscosity = scale*(1e6)*p[0]
viscosity_p = np.abs(viscosity-scale*(1e6)*p_upp[0])
viscosity_m = np.abs(viscosity-scale*(1e6)*p_low[0])
viscosity_pm = 0.5*(viscosity_p+viscosity_m)

print("-----------------------------------------------------")
print("Box volume =", V, "+/-", std_V, "nm^3")
print("-----------------------------------------------------")
print("eta =", viscosity, "+/-", viscosity_pm, "cP")
print("-----------------------------------------------------")

fig, (ax1,ax2) = plt.subplots(1,2)

ax1.plot(tint, avint[num_sets], 'k.')
ax1.plot(tfit, np.polyval(p,tfit), 'r--',
    linewidth=3)
for jj in range(num_sets) :
    ax1.plot(tint, avint[jj], 'b:',
    linewidth=0.75)
ax1.set_xlabel(r'$t$ [ps]', fontsize=20)
ax1.set_ylabel(r'$I^2$ [cP$^2$]', fontsize=20)
ax1.set_xlim([tint[0], tint[-1]])

ax2.plot(tint[1:], std_avint_p[1:], 'b-')
ax2.set_xlabel(r'$t$ [ps]', fontsize=20)
ax2.set_ylabel(r'$std(I^2)$ [$\%$]', fontsize=20)
ax2.set_xlim([tint[0], tint[-1]])

plt.show()