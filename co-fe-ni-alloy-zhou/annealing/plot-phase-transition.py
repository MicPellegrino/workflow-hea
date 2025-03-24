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

na = 108000

time, temperature_heating, n = read_xvg("temperature_heating.dat")
time, potential_heating, n = read_xvg("potential_heating.dat")
time, kinetic_heating, n = read_xvg("kinetic_heating.dat")

time, temperature_cooling, n = read_xvg("temperature_10.dat")
time, potential_cooling, n = read_xvg("potential_10.dat")
time, kinetic_cooling, n = read_xvg("kinetic_10.dat")

temperature_heating = np.array(temperature_heating[0])
potential_heating = np.array(potential_heating[0])/na
kinetic_heating = np.array(kinetic_heating[0])/na
total_heating = potential_heating+kinetic_heating

temperature_cooling = np.array(temperature_cooling[0])
potential_cooling = np.array(potential_cooling[0])/na
kinetic_cooling = np.array(kinetic_cooling[0])/na
total_cooling = potential_cooling+kinetic_cooling

# plt.plot(temperature_cooling,potential_cooling,'b.',label='Potential per-atom (cooling)')
# plt.plot(temperature_cooling,kinetic_cooling,'b.',label='Kinetic per-atom (cooling)')
plt.plot(temperature_cooling,kinetic_cooling+potential_cooling,'b.',label='Total per-atom (cooling)')
# plt.plot(temperature_heating,potential_heating,'ro',label='Potential per-atom (heating)')
# plt.plot(temperature_heating,kinetic_heating,'ro',label='Kinetic per-atom (heating)')
plt.plot(temperature_heating,kinetic_heating+potential_heating,'ro',label='Total per-atom (heating)')
plt.xlim([min(temperature_heating),max(temperature_heating)])
plt.legend(fontsize=25)
plt.tick_params(axis='both', labelsize=20)
plt.xlabel(r'$T$ [K]', fontsize=25)
plt.ylabel(r'$E$ [eV]', fontsize=25)
plt.show()