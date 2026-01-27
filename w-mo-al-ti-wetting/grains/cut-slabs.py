from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
import numpy as np
from copy import deepcopy

in_file_bulk = 'eq_5.34.data'
out_file_surface_1 = 'tests/eq_5d34_shift.data'
out_file_surface_2 = 'tests/eq_5d34_rotx.data'
out_file_surface_3 = 'tests/eq_5d34_roty.data'
out_file_droplet = 'tests/droplet.data'

pipeline_substrate = import_file(in_file_bulk)

data_substrate = pipeline_substrate.compute()
substrate_cell_0 = deepcopy(data_substrate.cell[...])
Lx = substrate_cell_0[0,0]
Ly = substrate_cell_0[1,1]
Lz = substrate_cell_0[2,2]
ox = substrate_cell_0[0,3]
oy = substrate_cell_0[1,3]
oz = substrate_cell_0[2,3]

shift1 = AffineTransformationModifier(
    transformation = [  [1, 0, 0, -ox],
                        [0, 1, 0, -oy],
                        [0, 0, 1, -oz]])
pipeline_substrate.modifiers.append(shift1)

data_substrate = pipeline_substrate.compute()
export_file(pipeline_substrate, out_file_surface_1, "lammps/data")

rot1 = AffineTransformationModifier(operate_on = {'particles'},
    transformation = [  [1, 0,  0, 0],
                        [0, 0, -1, Lz],
                        [0, 1,  0, 0]])
pipeline_substrate.modifiers.append(rot1)

data_substrate = pipeline_substrate.compute()
export_file(pipeline_substrate, out_file_surface_2, "lammps/data")

rot2 = AffineTransformationModifier(operate_on = {'particles'},
    transformation = [  [0, 0,  1,  0],
                        [0, 1,  0,  0],
                        [-1, 0, 0, Lx]])
pipeline_substrate.modifiers.append(rot2)

data_substrate = pipeline_substrate.compute()
export_file(pipeline_substrate, out_file_surface_3, "lammps/data")

n_slabs = 7
z_ref = 45.0
dz = Lz/n_slabs
for n in range(n_slabs) :
    data_substrate = pipeline_substrate.compute()
    data_substrate.apply(ExpressionSelectionModifier(expression = f'Position.Z<{n*dz} || Position.Z>={(n+1)*dz}'))
    data_substrate.apply(DeleteSelectedModifier(operate_on={'particles'}))
    z_com = np.mean(data_substrate.particles.positions[:,2])
    shift2 = AffineTransformationModifier(operate_on = {'particles'},
        transformation = [  [1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 1, z_ref-z_com]])
    data_substrate.apply(shift2)
    z_com = np.mean(data_substrate.particles.positions[:,2])
    export_file(data_substrate, f'tests/slab_{n}.data', "lammps/data")

droplet_file = '/home/michele/workflow-hea/ti-al-droplets/droplets/020/molten_AlTi_droplet_eq.data'
pipeline_droplet = import_file(droplet_file)

target1 = AffineTransformationModifier(operate_on = {'cell'},
    relative_mode = False,
    target_cell = [ [Lx, 0, 0, 0],
                    [0, Ly, 0, 0],
                    [0, 0, Lz, 0]])
pipeline_droplet.modifiers.append(target1)
data_droplet = pipeline_droplet.compute()

x_ref_d = 0.5*Lx
y_ref_d = 0.5*Ly
z_ref_d = 102.5
x_com_d = np.mean(data_droplet.particles.positions[:,0])
y_com_d = np.mean(data_droplet.particles.positions[:,1])
z_com_d = np.mean(data_droplet.particles.positions[:,2])
shift3 = AffineTransformationModifier(operate_on = {'particles'},
    transformation = [  [1, 0, 0, x_ref_d-x_com_d],
                        [0, 1, 0, y_ref_d-y_com_d],
                        [0, 0, 1, z_ref_d-z_com_d]])
data_substrate.apply(shift3)
pipeline_droplet.modifiers.append(shift3)
data_droplet = pipeline_droplet.compute()

export_file(pipeline_droplet, out_file_droplet, "lammps/data")