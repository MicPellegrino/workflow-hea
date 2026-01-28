from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
import numpy as np
from copy import deepcopy
import lammps

def align_3_axes(in_file_bulk, out_folder) :
    
    pipeline_bulk = import_file(in_file_bulk)
    data_bulk = pipeline_bulk.compute()
    cell_0 = deepcopy(data_bulk.cell[...])
    
    # Save dimensions and initial origin
    Lx = cell_0[0,0]
    Ly = cell_0[1,1]
    Lz = cell_0[2,2]
    ox = cell_0[0,3]
    oy = cell_0[1,3]
    oz = cell_0[2,3]

    # Shift so that the origin is now at (0,0,0)
    # For bloody reasons it's not possible to o the same with a target box!
    # Curse you OVITO!
    shift = AffineTransformationModifier(
        transformation = [[1, 0, 0, -ox],
                          [0, 1, 0, -oy],
                          [0, 0, 1, -oz]])
    pipeline_bulk.modifiers.append(shift)
    data_bulk = pipeline_bulk.compute()
    export_file(pipeline_bulk, out_folder+'/bulk_z.data', "lammps/data")

    # Rotate along the x-axis and output
    rotx = AffineTransformationModifier(operate_on = {'particles'},
        transformation = [[1, 0,  0, 0],
                          [0, 0, -1, Lz],
                          [0, 1,  0, 0]])
    pipeline_bulk.modifiers.append(rotx)
    data_bulk = pipeline_bulk.compute()
    export_file(pipeline_bulk, out_folder+'/bulk_x.data', "lammps/data")

    # Rotate along the y-axis and output
    roty = AffineTransformationModifier(operate_on = {'particles'},
        transformation = [[0, 0,  1,  0],
                          [0, 1,  0,  0],
                          [-1, 0, 0, Lx]])
    pipeline_bulk.modifiers.append(roty)
    data_bulk = pipeline_bulk.compute()
    export_file(pipeline_bulk, out_folder+'/bulk_y.data', "lammps/data")

def cut_slabs(in_file_bulk, out_folder, n_slabs, z_ref) :

    pipeline_bulk = import_file(in_file_bulk)
    data_bulk = pipeline_bulk.compute()
    cell_0 = deepcopy(data_bulk.cell[...])
    Lx = cell_0[0,0]
    Ly = cell_0[1,1]
    Lz = cell_0[2,2]

    # Cut slabs of width dz and shift their COM (or rather centre of geometry) 
    # to the reference z coordinate
    dz = Lz/n_slabs
    for n in range(n_slabs) :
        data_substrate = pipeline_bulk.compute()
        data_substrate.apply(ExpressionSelectionModifier(expression = f'Position.Z<{n*dz} || Position.Z>={(n+1)*dz}'))
        data_substrate.apply(DeleteSelectedModifier(operate_on={'particles'}))
        z_com = np.mean(data_substrate.particles.positions[:,2])
        shift = AffineTransformationModifier(operate_on = {'particles'},
            transformation = [[1, 0, 0,           0],
                              [0, 1, 0,           0],
                              [0, 0, 1, z_ref-z_com]])
        data_substrate.apply(shift)
        export_file(data_substrate, out_folder+f'/slab_{n}.data', "lammps/data")

def align_droplet(in_droplet_file, out_droplet_file, in_refconf_file, z_ref_d) :

    pipeline_ref = import_file(in_refconf_file)
    data_ref = pipeline_ref.compute()
    cell_0 = deepcopy(data_ref.cell[...])
    Lx = cell_0[0,0]
    Ly = cell_0[1,1]
    Lz = cell_0[2,2]

    pipeline_droplet = import_file(in_droplet_file)
    target = AffineTransformationModifier(operate_on = {'cell'},
        relative_mode = False,
        target_cell = [[Lx, 0, 0, 0],
                       [0, Ly, 0, 0],
                       [0, 0, Lz, 0]])
    pipeline_droplet.modifiers.append(target)
    data_droplet = pipeline_droplet.compute()

    x_ref_d = 0.5*Lx
    y_ref_d = 0.5*Ly
    x_com_d = np.mean(data_droplet.particles.positions[:,0])
    y_com_d = np.mean(data_droplet.particles.positions[:,1])
    z_com_d = np.mean(data_droplet.particles.positions[:,2])
    shift = AffineTransformationModifier(operate_on = {'particles'},
        transformation = [[1, 0, 0, x_ref_d-x_com_d],
                          [0, 1, 0, y_ref_d-y_com_d],
                          [0, 0, 1, z_ref_d-z_com_d]])
    pipeline_droplet.modifiers.append(shift)
    data_droplet = pipeline_droplet.compute()
    export_file(pipeline_droplet, out_droplet_file, "lammps/data")

def combine_with_lammps(in_file_droplet, in_file_substrate, out_file_combined, ff_file) :

    # Use a dummy LAMMPS code to merge the droplet and the substrate configuration.
    # This ensures that the output configuration is LAMMPS-readable and avoid 
    # messing up with the .data file format.
    lmp = lammps.lammps()
    read_data1 = in_file_droplet
    read_data2 = in_file_substrate
    write_data = out_file_combined
    lmp.commands_list([
        'units metal',
        'atom_style atomic',
        'pair_style eam/alloy',
        'neighbor 2.0 bin',
        'neigh_modify check yes',
        'boundary p p p'])
    c1 = "read_data \""+read_data1+"\" extra/atom/types 2"
    c2 = "read_data \""+read_data2+"\" add append offset 2 0 0 0 0 "
    lmp.command(c1)
    lmp.command(c2)
    lmp.commands_list([
        f'pair_coeff * * {ff_file} Al Ti W Mo',
        'group liquid type 1 2',
        'group solid type 3 4'])
    c3 = "write_data \""+write_data+"\" "
    lmp.command(c3)


if __name__ == "__main__" :
    
    print("### Aligning the bulk box along the three Cartesian axes ###")
    in_file_bulk = 'eq_5.34.data'
    out_folder = '5d34_confs'
    align_3_axes(in_file_bulk, out_folder)

    print("### Cutting solid slabs                                  ###")
    in_file_bulk = '5d34_confs/bulk_x.data'
    out_folder = '5d34_confs/slabs_x'
    cut_slabs(in_file_bulk, out_folder, n_slabs=7, z_ref=45.0)

    print("### Aligning droplet to match substrate cell dimensions  ###")
    in_droplet_file = '/home/michele/workflow-hea/ti-al-droplets/droplets/020/molten_AlTi_droplet_eq.data'
    out_droplet_file = '5d34_confs/droplet_020.data'
    in_refconf_file = '5d34_confs/bulk_x.data'
    align_droplet(in_droplet_file, out_droplet_file, in_refconf_file, z_ref_d=102.5)

    print("### Combine droplet and substrate configurations         ###")
    in_file_droplet = '5d34_confs/droplet_020.data'
    in_file_substrate = '5d34_confs/slabs_x/slab_1.data'
    out_file_combined = '5d34_confs/input_confs/combined_020_x1.data'
    ff_file = '/home/michele/workflow-hea/potential-files/eam/CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy'
    combine_with_lammps(in_file_droplet, in_file_substrate, out_file_combined, ff_file)