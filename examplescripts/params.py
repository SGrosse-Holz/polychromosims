import numpy as np
from mootils import paramtools

# Control center: a bunch of useful switches
steps_per_block = 1000 # TODO: clean this up
total_blocks = 10
save_to = "generic" # "new folder", "generic" or path to an existing folder

start_from = "custom" # cubic, conf, custom
add_spherical_confinement = True
add_lamina_attraction = True
add_pulling = False

# Saving data
base_folder = "/net/levsha/share/simongh/sims/data/locus_pulling_equilibration"
folder_flag = "equil" # arbitrary string to add to the folder name
blocks_per_file = 100

# Computational parameters
# For reference, see polychrom/simulation.py
platform = "CUDA"
GPU = "0" # Note: might be overridden by command line to sim.py
integrator = "variableLangevin"
collision_rate = 0.01
error_tol = 0.001   # for variableLangevin
timestep = 20      # for Brownian
verbose = False
max_Ek = 10
PBCbox = False # [x, y, z] or False or "density" to use sphconf_density

# Major simulation parameters
N = 60000 # Note: this will be overwritten if loading a conformation

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

# Lamina attraction
lam_r = "sphconf_density" # "sphconf_density" or number. The former will infer from spherical confinement.
lam_width = 2.0
lam_depth = 3.0
lam_particles = paramtools.lamina_particles()

# Polymer chain forcekit
chains = [(0, None, False)] # Likely to be overwritten
chain_bondlength = 1.0
chain_bondwiggledist = 0.05
chain_anglek = 1.5 # From example.py: 1.5 is reasonable, 4 gives l_p = 4mon, 8 is very stiff
chain_nonbonded = "hetero" # "hetero", "simple"
chain_nb_extrahard = "all" # "all", list of indices, or a fraction between 0 and 1, to be assigned at random
chain_nb_repEnergy = 3.0
chain_nb_repEnergy_extrahard = 20.0
chain_nb_repRadius = 1.0
chain_nb_attrEnergy = 0.0
chain_nb_attrRadius = 2.0

hetero_energy = 1.0
hetero_matrix = np.array([[0.0, 0.0, 0.0],
                          [0.0, 0.4, 0.6],
                          [0.0, 0.0, 1.0]])
hetero_monTypes = paramtools.gen_compartments_Ed(pulled_locus_is_C = True)
hetero_overwrite_conf = False # If there are compartment IDs in the given starting
                        # conformation, should they be overwritten with this?

# Locus pulling
pulled_locus = paramtools.pulled_locus_Ed()
pull_direction = "center" # [x, y, z] or "center": pull towards center of mass
pull_force = 1.0

# Energy minimization
do_energy_minimization = True # Automatically switched off when starting from conf
emin_tolerance=0.3
emin_maxIter = 0
emin_randomOffset = chain_bondwiggledist / 2
