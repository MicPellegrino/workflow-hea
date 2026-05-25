import numpy as np 
import matplotlib.pyplot as plt 
from read_profiles import read_scalar_field, merge_fields

LABEL_FS = 25
LEGEND_FS = 20
TICKS_FS = 20
MARKER_SIZE = 12

FRAME_TO_PLOT = -1

a = read_scalar_field("slab_position_check_sol.txt", "bin_loc")
b = read_scalar_field("temp_chunk_bias_sol.txt", "temperature")
data_sol = merge_fields(a, b)

a = read_scalar_field("slab_position_check_liq.txt", "bin_loc")
b = read_scalar_field("temp_chunk_bias_liq.txt", "temperature")
data_liq = merge_fields(a, b)

frames = data_sol.keys()

assert frames==data_liq.keys(), "Solid and liquid temperature profiles computed on different frames!"

frames = list(frames)

bins_liq = data_liq[frames[FRAME_TO_PLOT]]['bin_loc']
temp_liq = data_liq[frames[FRAME_TO_PLOT]]['temperature']
bins_sol = data_sol[frames[FRAME_TO_PLOT]]['bin_loc']
temp_sol = data_sol[frames[FRAME_TO_PLOT]]['temperature']
plt.plot(bins_liq,temp_liq,'ro',ms=MARKER_SIZE,label="liquid")
plt.plot(bins_sol,temp_sol,'bx',ms=MARKER_SIZE,label="solid")
plt.legend(fontsize=LEGEND_FS)
plt.xlabel('z [A]', fontsize=LABEL_FS)
plt.ylabel('T [K]', fontsize=LABEL_FS)
plt.xticks(fontsize=TICKS_FS)
plt.yticks(fontsize=TICKS_FS)
plt.show()

n_frames = len(frames)

for i in range(n_frames) :
    bins_liq = data_liq[frames[i]]['bin_loc']
    temp_liq = data_liq[frames[i]]['temperature']
    bins_sol = data_sol[frames[i]]['bin_loc']
    temp_sol = data_sol[frames[i]]['temperature']
    alpha = ((i+1)/n_frames)**2
    if i == (n_frames-1) :
        plt.plot(bins_liq,temp_liq,'ro',ms=MARKER_SIZE,alpha=alpha,label="liquid")
        plt.plot(bins_sol,temp_sol,'bx',ms=MARKER_SIZE,alpha=alpha,label="solid")
    else :
        plt.plot(bins_liq,temp_liq,'ro',ms=MARKER_SIZE,alpha=alpha)
        plt.plot(bins_sol,temp_sol,'bx',ms=MARKER_SIZE,alpha=alpha)
plt.legend(fontsize=LEGEND_FS)
plt.xlabel('z [A]', fontsize=LABEL_FS)
plt.ylabel('T [K]', fontsize=LABEL_FS)
plt.xticks(fontsize=TICKS_FS)
plt.yticks(fontsize=TICKS_FS)
plt.show()