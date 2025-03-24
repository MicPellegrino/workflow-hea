import numpy as np 
import matplotlib.pyplot as plt 

t_cool = [1,4,5,8,10]

phase_frac = {
    'FCC': [],
    'BCC': [],
    'HCP': [],
    'AMO': [],
}

a1 = 0.333
a2 = 0.334
a3 = 0.333

labels = []

for t in t_cool :

    data = np.load('annealed_system_'+str(t)+'.npz')
    f_fcc = data['f_fcc_vec']
    f_bcc = data['f_bcc_vec']
    f_hcp = data['f_hcp_vec']
    f_amo = data['f_amo_vec']
    
    labels.append(str(t))

    phase_frac['FCC'].append(f_fcc[-1])
    phase_frac['BCC'].append(f_bcc[-1])
    phase_frac['HCP'].append(f_hcp[-1])
    phase_frac['AMO'].append(f_amo[-1])

    # plt.hist(t,f_amo[-1],'ro')
    # plt.hist(t,f_hcp[-1],'gd')
    # plt.hist(t,f_fcc[-1],'ks')
    # plt.hist(t,f_bcc[-1],'bx')

width=0.5
bottom = np.zeros(len(labels))
fig, ax = plt.subplots()

for phase, frac in phase_frac.items():
    p = ax.bar(labels, frac, width, label=phase, bottom=bottom)
    bottom += np.array(frac)
    ax.bar_label(p, label_type='center')

ax.set_title('CoFeNi bulk phases')
ax.set_xlabel('Cooling time [ns]')
ax.set_ylabel('Fraction')
ax.legend()
plt.show()

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
file_path = 'temperature_5.dat'
columns = read_dat_file(file_path)

plt.plot(dt*np.array(columns[0]),columns[1],'k-')
plt.xlabel('t [ps]')
plt.ylabel('T [K]')
plt.show()