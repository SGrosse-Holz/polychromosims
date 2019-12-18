import numpy as np
from mootils import paramtools

# Control center: a bunch of useful switches
steps_per_block = 1000
total_blocks = 10
save_to = "generic" # "new folder", "generic" or path to an existing folder

start_from = "cubic" # cubic, conf, custom
add_spherical_confinement = True
add_lamina_attraction = False
add_pulling = False
add_crosslinks = False
has_extrusion = False

# Saving data
base_folder = "/net/levsha/share/simongh/sims/data"
folder_flag = "example" # arbitrary string to add to the folder name
blocks_per_file = 100

# Computational parameters
# For reference, see polychrom/simulation.py
platform = "CUDA"
GPU = "0" # likely to be overridden
integrator = "variableLangevin"
collision_rate = 0.03
error_tol = 0.01   # for variableLangevin
timestep = 20      # for Brownian
verbose = False
max_Ek = 10
PBCbox = False # [x, y, z] or False or "density" to use sphconf_density

# Major simulation parameters
N = 10000 # Note: this will be overwritten if loading a conformation

# Starting conformation
start_cubic_box = None # None if should be computed from N

### FORCES ###
# Spherical confinement
sphconf_density = 0.35
sphconf_k = 6

# Polymer chain forcekit
chains = [(0, None, False)] # Likely to be overwritten
chain_bondlength = 1.0
chain_bondwiggledist = 0.05
chain_anglek = 1.5 # From example.py: 1.5 is reasonable, 4 gives l_p = 4mon, 8 is very stiff

chain_nonbonded = "polynomial" # "hetero", "simple", "polynomial"
chain_nb_repEnergy = 3.0

chain_nb_repRadius = 1.0
chain_nb_attrEnergy = 0.0
chain_nb_attrRadius = 2.0

# Energy minimization
do_energy_minimization = True # Automatically switched off when starting from conf
emin_tolerance=0.3
emin_maxIter = 0
emin_randomOffset = chain_bondwiggledist / 2
