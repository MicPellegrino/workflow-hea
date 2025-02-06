from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

import numpy as np

from copy import deepcopy

# Set shift_z to a different number and/or
# define a shift (new mod) for the substrate too
scale_z = 0.75
shift_z = 23

pipeline_droplet = import_file("droplet_Al.data")
pipeline_substrate = import_file("annealed_surface_WMo.data")

data_droplet = pipeline_droplet.compute()
droplet_cell_0 = deepcopy(data_droplet.cell[...])

data_substrate = pipeline_substrate.compute()
substrate_cell_0 = deepcopy(data_substrate.cell[...])

Lx = substrate_cell_0[0,0]
Ly = substrate_cell_0[1,1]
Lz = scale_z*droplet_cell_0[2,2]

ax = Lx/droplet_cell_0[0,0]
ay = Ly/droplet_cell_0[1,1]

dx = 0.5*(Lx-droplet_cell_0[0,0])-droplet_cell_0[0,3]
dy = 0.5*(Ly-droplet_cell_0[1,1])-droplet_cell_0[1,3]
dz = -droplet_cell_0[2,3]

print("Initial droplet box:")
print(droplet_cell_0)

cell_matrix = np.array(droplet_cell_0)
cell_matrix[:,0] = np.array([Lx,0.0,0.0])
cell_matrix[:,1] = np.array([0.0,Ly,0.0])
cell_matrix[:,2] = np.array([0.0,0.0,Lz])
cell_matrix[:,3] = np.array([0.0,0.0,0.0])

print(cell_matrix)

mod11 = AffineTransformationModifier(operate_on = {'cell'},
			transformation = [[1, 0, 0, -droplet_cell_0[0,3]],
                              [0, 1, 0, -droplet_cell_0[1,3]],
                              [0, 0, 1, -droplet_cell_0[2,3]]])

mod12 = AffineTransformationModifier(operate_on = {'cell'},
			transformation = [[ax, 0, 0, 0],
                              [0, ay, 0, 0],
                              [0, 0, scale_z, 0]])

# mod1 = AffineTransformationModifier(reduced_coords=False, target_cell=cell_matrix)

pipeline_droplet.modifiers.append(mod11)
data_droplet = pipeline_droplet.compute()
pipeline_droplet.modifiers.append(mod12)
data_droplet = pipeline_droplet.compute()

mod2 = AffineTransformationModifier(operate_on = {'particles'},
                              	transformation = [[1, 0, 0, dy],
                                                  [0, 1, 0, dx],
                                                  [0, 0, 1, dz]])

pipeline_droplet.modifiers.append(mod2)
data_droplet = pipeline_droplet.compute()

mod3 = AffineTransformationModifier(operate_on = {'particles'},
                                transformation = [[1, 0, 0, 0],
                                                  [0, 1, 0, 0],
                                                  [0, 0, 1, -shift_z]])
pipeline_droplet.modifiers.append(mod3)
data_droplet = pipeline_droplet.compute()

droplet_cell_1 = deepcopy(data_droplet.cell[...])

print("Final droplet box:")
print(droplet_cell_1)
export_file(pipeline_droplet, "droplet_Al-transformed.data", "lammps/data")

print("Initial substrate box:")
print(substrate_cell_0)

mod4 = AffineTransformationModifier(transformation = [[1, 0, 0, -substrate_cell_0[0,3]],
                                                  [0, 1, 0, -substrate_cell_0[1,3]],
                                                  [0, 0, 1, -substrate_cell_0[2,3]]])
pipeline_substrate.modifiers.append(mod4)
data_substrate = pipeline_substrate.compute()

substrate_cell_1 = deepcopy(data_substrate.cell[...])

print("Final substrate box:")
print(substrate_cell_1)

export_file(pipeline_substrate, "annealed_surface_WMo-transformed.data", "lammps/data")
