from ovito.io import *
from ovito.modifiers import *

import numpy as np

pipeline_wetting = import_file("wetting.dump")
pipeline_wetting.modifiers.append(CommonNeighborAnalysisModifier())

for data in pipeline_wetting.frames:
    
    # TEST #
    # n_bcc = data.attributes['CommonNeighborAnalysis.counts.BCC']
    # print(n_bcc / data.particles.count)
    
    # Example: average z corrdinate
    coord_liq = data.particles.positions[data.particles['v_dummymol'] == 1]
    coord_liq = coord_liq[...]
    zmean_liq = np.mean(coord_liq[:,2])
    print(zmean_liq)
    coord_sub = data.particles.positions[data.particles['v_dummymol'] == 2]
    coord_sub = coord_sub[...]
    zmean_sub = np.mean(coord_sub[:,2])
    print(zmean_sub)
    
