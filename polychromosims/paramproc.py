import os,sys,shutil
import importlib
import datetime

import numpy as np

from polychrom.hdf5_format import load_hdf5_file

from mootils.save_module_to_script import mod2py
import polychromosims.globalvars as globalvars

def update_from_cmd(params):
    # First things first: parse cmd arguments
    sys.stderr.write(str(sys.argv)+"\n")
    for arg in sys.argv:
        if arg[:7] == "params.":
            sys.stderr.write("setting "+arg+"\n")
            exec(arg)

def proc(params):
    # ------------------- Resolving special keywords in the parameters ----------------------------
    # NOTE: this should not set anything unconditionally, such that it does not
    # unexpectedly overwrite anything. It is just supposed to resolve keywords
    # to actual parameter values.

    # Set conformation parameters
    # This might overwrite certain parameters (almost certainly N), so it should be
    # done before other processing.
    if params.start_from == "conf":
        params.do_energy_minimization = False
        conf = load_hdf5_file(params.start_conf_name)
        try:
            params.N = len(conf['pos'])
            chains = conf['chains']
            params.chains = [(chains[i, 0], chains[i, 1], bool(chains[i, 2])) for i in range(chains.shape[0])]
            if params.chain_nonbonded == "hetero" and not params.hetero_overwrite_conf:
                try:
                    params.hetero_monTypes = conf['hetero_monTypes']
                except KeyError: # Backwards compatibility
                    params.hetero_monTypes = conf['c5_compartmentIDs']
        except KeyError as key:
            print("WARNING: Did not find "+str(key)+" in file "+params.start_conf_name+". Using default.")
    elif params.start_from == "cubic" and params.start_cubic_box == None:
        params.start_cubic_box = 5*int(params.N**(1/3))
    elif params.start_from == "custom":
        spec = importlib.util.spec_from_file_location("startconf", params.start_custom_script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        for var in  [var for var in dir(mod) if not var[0:2] == "__"]:
            setattr(params, var, getattr(mod, var))

    if type(params.hetero_monTypes) == np.ndarray:
        params.hetero_monTypes = params.hetero_monTypes.astype(int)
    
    # If desired, set PBC box according to given density
    if params.PBCbox == "density":
        pbcdim = (params.N * params.sphconf_density)**(1/3)
        params.PBCbox = [pbcdim, pbcdim, pbcdim]
    
    # Convert chain_nb_extrahard to an actual list
    if isinstance(params.chain_nb_extrahard, str) and params.chain_nb_extrahard == "all":
        params.chain_nb_extrahard = np.arange(params.N)
    elif isinstance(params.chain_nb_extrahard, (int, float)):
        params.chain_nb_extrahard = np.random.choice(params.N, int(params.chain_nb_extrahard*params.N), replace=False)
    # else: chain_nb_extrahard is a list of indices already
    
    # Lamina attraction radius
    if params.lam_r == "sphconf_density":
        params.lam_r = (3 * params.N / (4 * 3.141592 * params.sphconf_density)) ** (1 / 3.) # from polychrom.forces.spherical_confinement
    
    # Generate the folder name and create it
    if params.save_to == "new folder":
        now = datetime.datetime.now()
        folder = os.path.join(params.base_folder,
                "{0:04}-{1:02}-{2:02}_{3:02}:{4:02}".format(now.year,
                                                            now.month,
                                                            now.day,
                                                            now.hour,
                                                            now.minute)
                )
        folder += "_" + params.folder_flag
        rep = 0
        repfolder = lambda : folder+"_"+str(rep)
        while globalvars.params_allow_execution:
            try:
                os.makedirs(repfolder())
            except FileExistsError:
                rep += 1
                continue
            break
        params.folder = repfolder()
    elif params.save_to == "generic":
        params.folder = "/net/levsha/share/simongh/sims/data/generic_trajectory"
        if globalvars.params_allow_execution:
            shutil.rmtree(params.folder)
            os.makedirs(params.folder)
    else:
        params.folder = params.save_to
    # Make sure that upon a repeated execution we do not unintentionally change the
    # folder name
    params.save_to = params.folder
    
def write_processed(params):
    bare_filename = os.path.join(params.folder, 'processed_params')
    rep = 0
    repfile = lambda : bare_filename+'_'+str(rep)+'.py'
    while True:
        try:
            mod2py(params, repfile())
        except FileExistsError:
            rep += 1
            continue
        break

def start_editing(params):
    # Save the stuff that should not be reloaded
    folder = params.folder
    extrahard = params.chain_nb_extrahard

    params_modspec = globalvars.params_modspec

    # Reload
    importlib.reload(globalvars)
    params_modspec.loader.exec_module(params)
    globalvars.params_modspec = params_modspec

    # Set saved stuff
    params.save_to = folder
    params.chain_nb_extrahard = extrahard

    # Get command line values
    update_from_cmd(params)
    sys.stderr.write("After start_editing: GPU = "+params.GPU+"\n")

def end_editing(params):
    proc(params)
    sys.stderr.write("After editing params: GPU = "+params.GPU+"\n")
