import os,sys,shutil
import importlib
import warnings
import datetime

import polychromosims.paramproc
import polychromosims.globalvars as globalvars

def init(sysargv=sys.argv):
    """
    Initialize the polychromosims wrapper. This function should be called
    before importing polychromosims.sim, because it loads the params module.

    Parameters
    ----------
    sysargv : list
        a list of parameters, default is sys.argv

    Returns
    -------
    the loaded params module
    """
    # The paramscript is given as commandline argument, so we have to load it dynamically.
    # NOTE: the GPU=... argument was replaced by giving params.GPU = ...
    myparamsfile = (["params.py"] + [arg[12:] for arg in sysargv if 'paramscript=' in arg])[-1]
    params_modspec = importlib.util.spec_from_file_location("params", myparamsfile)
    params = importlib.util.module_from_spec(params_modspec)
    params_modspec.loader.exec_module(params)
    sys.modules[params.__name__] = params
    globalvars.params_modspec = params_modspec
    
    polychromosims.paramproc.update_from_cmd(params)
    globalvars.params_allow_execution = True
    polychromosims.paramproc.proc(params)
    globalvars.params_module = params
    import polychromosims.sim as sim

    # Save relevant files with the data
    shutil.copyfile(myparamsfile, os.path.join(params.folder, 'params.py'))
    shutil.copyfile(os.path.abspath(sysargv[0]), os.path.join(params.folder, 'simscript.py'))

    return params, sim

def finalize():
    print("saved everything to "+globalvars.params_module.folder)
    print("It is now "+str(datetime.datetime.now()))
