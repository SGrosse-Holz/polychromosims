import os,sys
import numpy as np

from polychrom.hdf5_format import HDF5Reporter

import polychromosims
from polychromosims.paramproc import start_editing, end_editing, write_processed

# Set up
params, sim = polychromosims.init(sys.argv)
reporter = HDF5Reporter(folder=params.folder,
        max_data_length=params.blocks_per_file, overwrite=True)

# Equilibration (pretty pseudo here)
start_editing(params)

params.integrator = "variableLangevin"
params.collision_rate = 0.01
params.steps_per_block = 1000
params.total_blocks = 2

end_editing(params)
sim.fullSim(reporter)

# Run
start_editing(params)
params.start_from = "conf"
params.start_conf_name = os.path.join(params.folder, "final_conformation_0.h5")

end_editing(params)
write_processed(params)
sim.fullSim(reporter)

# Clean-up
polychromosims.finalize()
exit()
