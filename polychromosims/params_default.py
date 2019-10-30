# This file serves a double purpose:
#  - it's an example parameter script, showing all the variables that can be
#    set
#  - it has the default values, such that we don't have to take care to include
#    everything in any specific script
# 
# As such, I'll try to always keep it up to date with developments of sim.py
# and paramproc.py

import numpy as np
from mootils import paramtools

# Control center: a bunch of useful switches
steps_per_block = 5000
total_blocks = 10
save_to = "generic" # "new folder", "generic" or path to an existing folder

start_from = "cubic" # cubic, conf, custom
add_spherical_confinement = False
add_lamina_attraction = False
add_pulling = False
add_crosslinks = False
has_extrusion = False

# Saving data
base_folder = "/net/levsha/share/simongh/sims/data"
folder_flag = "!defaults!" # arbitrary string to add to the folder name
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
start_conf_name = "/net/levsha/share/simongh/sims/startconf/locus_pulling/repeat_0.h5"
start_custom_script = "/net/levsha/share/simongh/sims/startconf/cubicComps.py" # a module defining
        # a) a function gen_start_conf that actually gives the conformation
        # b) all the meta params it needs, usually N, chains, hetero_monTypes
        #    and the like

### FORCES ###
# Spherical confinement
sphconf_density = 0.35
sphconf_k = 6

# Locus pulling
pulled_locus = []
pull_direction = "center" # [x, y, z] or "center": pull towards center of mass
pull_force = 0.2

# Crosslinks
cl_bonds = []

# Lamina attraction
lam_r = "sphconf_density" # "sphconf_density" or number. The former will infer from spherical confinement.
lam_width = 2.0
lam_depth = 3.0
lam_particles = []

# Polymer chain forcekit
chains = [(0, None, False)] # Likely to be overwritten
chain_bondlength = 1.0
chain_bondwiggledist = 0.05
chain_anglek = 1.5 # From example.py: 1.5 is reasonable, 4 gives l_p = 4mon, 8 is very stiff

chain_nonbonded = "polynomial" # "hetero", "simple", "polynomial"
chain_nb_repEnergy = 3.0

chain_nb_extrahard = [] # "all", list of indices, or a fraction between 0 and 1, to be assigned at random
chain_nb_repEnergy_extrahard = 20.0
chain_nb_repRadius = 1.0
chain_nb_attrEnergy = 0.0
chain_nb_attrRadius = 2.0

hetero_energy = 1.0
hetero_matrix = np.array([[0.0, 0.0, 0.0],
                          [0.0, 0.4, 0.0],
                          [0.0, 0.0, 1.0]])
hetero_monTypes = [] # This will throw an error if chain_nonbonded = "hetero"
hetero_overwrite_conf = False # If there are compartment IDs in the given starting
                        # conformation, should they be overwritten with this?

### EXTRUSION
# LEF parameters
extrusion_separation = 1000
extrusion_processivity = 100
extrusion_processivity_stalled = 100
extrusion_MDstepsPerStep = 1000

# CTCF parameters
extrusion_CTCFdict = {'captureLeft' : {}, 'captureRight' : {}, 'releaseLeft' : {}, 'releaseRight' : {}}
extrusion_p_capture = 1
extrusion_p_releasePerStep = 0

# Computational stuff
extrusion_stepsPerRestart = 100
extrusion_smc_bondwiggledist = 0.2
extrusion_smc_bondDist = 0.5

# Energy minimization
do_energy_minimization = True # Automatically switched off when starting from conf
emin_tolerance=0.3
emin_maxIter = 0
emin_randomOffset = chain_bondwiggledist / 2
