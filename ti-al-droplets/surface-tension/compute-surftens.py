import matplotlib.pyplot as plt 
import numpy as np

LMP_TO_SI = 1e-5

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

folder_name = '100'
dt = 0.001
t0 = 50

data = read_two_column_file(folder_name+'/surftens.dat')
data = np.array(data)
t = dt*data[:,0]
g = data[:,1]

n0 = np.argmin(np.abs(t-t0))

# In SI units
avg_surftens = LMP_TO_SI*np.mean(g[n0:])
std_surftens = LMP_TO_SI*np.std(g[n0:])

print("#####",folder_name,"#####")
print("<surf_tens> =", avg_surftens, "[Pa*m]")
print("std(surf_tens) =", std_surftens, "[Pa*m]")

"""
    Plotting test cases
"""
# labels = ['cpu','opencl_sp','opencl_dp','cuda_sp','cuda_dp']
# for l in labels :
#     data = read_two_column_file(folder_name+'_'+l+'/surftens.dat')
#     data = np.array(data)
#     t = dt*data[:,0]
#     g = data[:,1]
#     plt.plot(t,g,linewidth=3,label=l)
# plt.xlim([0,200])
# plt.legend(fontsize=25)
# plt.ylabel('surf. tens. [bar*A]',fontsize=25)
# plt.xlabel('time [ps]',fontsize=25)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.show()
