# This follows closely Max's example.py, but includes modifications and additions
# Note: this module should not be imported directly. Use polychromosims.initialize.init

from __future__ import absolute_import, division, print_function, unicode_literals
import os,sys
import shutil
import warnings
import importlib
import polychrom
from polychrom import (simulation, starting_conformations,
                       forces, forcekits)
import simtk.openmm as openmm
from polychrom.hdf5_format import HDF5Reporter, load_hdf5_file, load_URI

import numpy as np

import polychromosims.globalvars as globalvars
try:
    params = globalvars.params_module
except:
    raise RuntimeError("sim.py does not run alone! It needs a loaded"
                       "parameter module in globalvars.params_module.")

def createSim(reporters):
    """"""
    if not isinstance(reporters, list):
        reporters = [reporters]

    #Simulation object has many parameters that should be described in polychrom/simulation.py file 
    return simulation.Simulation(
            platform=params.platform,
            GPU = params.GPU,
            integrator=params.integrator,
            error_tol=params.error_tol,
            timestep=params.timestep,
            collision_rate=params.collision_rate,
            verbose=params.verbose,
            max_Ek=params.max_Ek,
            PBCbox=params.PBCbox,
            N = params.N,
            reporters=reporters) 
    
def setStartingConf(sim):
    """"""
    if params.start_from == "conf":
        polymer = load_hdf5_file(params.start_conf_name)['pos']
    elif params.start_from == "cubic":
        polymer = starting_conformations.grow_cubic(params.N, params.start_cubic_box)
    elif params.start_from == "custom":
        polymer = params.gen_start_conf()
    else:
        raise ValueError(params.errorString("start_from", params.start_from))
    
    sim.set_data(polymer, center=True)  # Do we want centering as an option?
    return sim

def addForces(sim):
    """"""
    # Confinement
    if params.add_spherical_confinement:
        sim.add_force( forces.spherical_confinement(sim,
            density=params.sphconf_density, k=params.sphconf_k))
    
    # Nonbonded force, possibly with compartments
    if params.chain_nonbonded == "hetero":
        nb_force = forces.heteropolymer_SSW
        nb_kwargs = {
                'interactionMatrix':params.hetero_matrix,
                'monomerTypes':params.hetero_monTypes,
                'extraHardParticlesIdxs':params.chain_nb_extrahard,
                'repulsionEnergy':params.chain_nb_repEnergy,
                'repulsionRadius':params.chain_nb_repRadius,
                'attractionEnergy':params.chain_nb_attrEnergy,
                'attractionRadius':params.chain_nb_attrRadius,
                'selectiveRepulsionEnergy':params.chain_nb_repEnergy_extrahard,
                'selectiveAttractionEnergy':params.hetero_energy,
                'name':'heteropolymer_SSW'}
    elif params.chain_nonbonded == "simple":
        nb_force = forces.smooth_square_well
        nb_kwargs = {
                'repulsionEnergy':params.chain_nb_repEnergy,
                'repulsionRadius':params.chain_nb_repRadius,
                'attractionEnergy':params.chain_nb_attrEnergy,
                'attractionRadius':params.chain_nb_attrRadius,
                'name':'smooth_square_well'}
    else:
        raise ValueError(params.errorString("chain_nonbonded", params.chain_nonbonded))

    sim.add_force(
        forcekits.polymer_chains(
            sim,
            chains=params.chains,
            bond_force_kwargs={
                'bondLength':params.chain_bondlength,
                'bondWiggleDistance':params.chain_bondwiggledist,
             },

            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                'k':params.chain_anglek, 
            },

            nonbonded_force_func=nb_force,
            nonbonded_force_kwargs=nb_kwargs,

            except_bonds=True,
        )
    )

    # Lamina attraction
    if params.add_lamina_attraction:
        sim.add_force(
            forces.spherical_well(
                sim,
                params.lam_particles,
                params.lam_r,
                width=params.lam_width,
                depth=params.lam_depth
            )
        )

    # Locus pulling
    if params.add_pulling:
        if not isinstance(params.pull_direction, list) and params.pull_direction == "center":
            com = np.mean(sim.get_data(), axis=0)
            com_locus = np.mean(sim.get_data()[params.pulled_locus], axis=0)
            params.pull_direction = com - com_locus
        pull_force = params.pull_force*params.pull_direction/np.linalg.norm(params.pull_direction)

        sim.add_force(
            forces.pull_force(
                sim,
                params.pulled_locus,
                [pull_force]
            )
        )

    return sim

def runSim(sim):
    if params.do_energy_minimization:
        sim.local_energy_minimization(tolerance=params.emin_tolerance,
                                      maxIterations=params.emin_maxIter,
                                      random_offset=params.emin_randomOffset)

    for _ in range(params.total_blocks):
        sim.do_block(params.steps_per_block)
    sim.print_stats()

    return sim

def saveFinal(sim, reporter):
    reportDict = {"pos":sim.get_data()}

    # Add chain information
    npchains = np.array(params.chains)
    npchains[np.equal(npchains, None)] = params.N
    reportDict['chains'] = npchains.astype(int)

    # Possibly add information about compartments
    if params.chain_nonbonded == "hetero":
        reportDict['hetero_monTypes'] = params.hetero_monTypes

    reporter.report("final_conformation", reportDict)

def fullSim(reporter):
    """
    A wrapper for everything happening here, assuming that all of this is
    basically always the same.
    """
    simobj = createSim(reporter)
    simobj = setStartingConf(simobj)
    simobj = addForces(simobj)
    simobj = runSim(simobj)
    saveFinal(simobj, reporter)
    reporter.dump_data()
