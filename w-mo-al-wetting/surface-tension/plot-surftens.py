import matplotlib.pyplot as plt 
import numpy as np

def read_two_column_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            columns = line.strip().split()
            if len(columns) == 2:
                data.append((float(columns[0]), float(columns[1])))
    return data

folder_name = 'slab_8_2'
labels = ['cpu','opencl_sp','opencl_dp','cuda_sp','cuda_dp']

dt = 0.001

for l in labels :
    data = read_two_column_file(folder_name+'_'+l+'/surftens.dat')
    data = np.array(data)
    t = dt*data[:,0]
    g = data[:,1]
    plt.plot(t,g,linewidth=3,label=l)
plt.xlim([0,200])
plt.legend(fontsize=25)
plt.ylabel('surf. tens. [bar*A]',fontsize=25)
plt.xlabel('time [ps]',fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()