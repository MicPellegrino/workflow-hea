import matplotlib.pyplot as plt
import numpy as np

dt = 0.001

def read_dat_file(file_path):
    columns = []
    with open(file_path, 'r') as file:
        lines = [line for line in file.readlines() if not line.startswith('#')]
        for line in lines:
            values = line.split()
            for i, value in enumerate(values):
                if len(columns) <= i:
                    columns.append([])
                columns[i].append(float(value))

    return columns

# Example usage:
file_path = 'temperature.dat'
columns = read_dat_file(file_path)

plt.plot(dt*np.array(columns[0]),columns[1],'k-')
plt.xlabel('t [ps]')
plt.ylabel('T [K]')
plt.show()