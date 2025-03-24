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

time, temperature, n = read_xvg("temperature_heating.dat")
time, potential, n = read_xvg("potential_heating.dat")
time, kinetic, n = read_xvg("kinetic_heating.dat")

temperature = np.array(temperature[0])
potential = np.array(potential[0])
kinetic = np.array(kinetic[0])
total = potential+kinetic

plt.plot(temperature,potential,'b.',label='Potential')
plt.plot(temperature,kinetic,'r.',label='Kinetic')
plt.plot(temperature,kinetic+potential,'k.',label='Total')
plt.legend()
plt.xlim([min(temperature),max(temperature)])
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$E$ [eV]')
plt.show()