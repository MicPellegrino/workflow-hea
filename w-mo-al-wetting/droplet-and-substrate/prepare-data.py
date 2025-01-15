from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

import numpy as np

from copy import deepcopy

pipeline_droplet = import_file("droplet_Al.data")
pipeline_substrate = import_file("annealed_surface_WMo.data")

data_droplet = pipeline_droplet.compute()
droplet_cell_0 = deepcopy(data_droplet.cell[...])

data_substrate = pipeline_substrate.compute()
substrate_cell_0 = deepcopy(data_substrate.cell[...])

substrate_cell_0[2,2] = droplet_cell_0[2,2]
print(droplet_cell_0)

# mod = AffineTransformationModifier(target_cell = substrate_cell_0)
mod1 = AffineTransformationModifier(transformation = [  [1, 0, 0, -droplet_cell_0[0,3]],
                                                        [0, 1, 0, -droplet_cell_0[1,3]],
                                                        [0, 0, 1, -droplet_cell_0[2,3]]])
pipeline_droplet.modifiers.append(mod1)
data_droplet = pipeline_droplet.compute()

mod2 = AffineTransformationModifier(operate_on = {'cell'},
                                    transformation = [  [1, 0, 0, 0],
                                                        [0, 1, 0, 0],
                                                        [0, 0, 0.75, 0]])
pipeline_droplet.modifiers.append(mod2)
data_droplet = pipeline_droplet.compute()

dz = 27.5
mod2 = AffineTransformationModifier(operate_on = {'particles'},
                                transformation = [  [1, 0, 0, 0],
                                                    [0, 1, 0, 0],
                                                    [0, 0, 1, -dz]])
pipeline_droplet.modifiers.append(mod2)
data_droplet = pipeline_droplet.compute()

droplet_cell_1 = deepcopy(data_droplet.cell[...])
print(droplet_cell_1)

export_file(pipeline_droplet, "droplet_Al-transformed.data", "lammps/data")