# polychromosims
My attempt at a useful wrapper for polychrom (which is an OpenMM wrapper...
where does this end?)

## Quick description

This is a wrapper for [polychrom](https://github.com/mirnylab/polychrom), the
Mirnylab's polymer simulation library. I use this to implement mechanisms that
are specific to projects but are still valuable to keep around. Accordingly,
this should be regarded more as a personal library; if you want to use it, you
will probably add your own features.

The basic idea of this package is to have one script ("simscript") that
specifies order of operations for a simulation, and another one ("paramscript")
that serves as a library/control center for the simulation parameters.

## Some details on where to find what

Implementation of features happens in `sim.py`, so this is the place to look
when wondering what a specific parameter actually does.

`params_default.py` can be used as a template for the parameter file. It
contains all the parameters that are used throughout the package. When modifying
anything, you can use `checks/params_default_completeness.py` to check that all
the parameters you use indeed do have default values.

## Examples

The `examples` folder contains everything that's needed to run a short
simulation with this wrapper. It even has the call to python wrapped into a
shell script, mostly for documentation purposes. So you should be able to run
this example by simply executing
```sh
$ ./run.sh
```

Note that the `params.py` in this folder is shortened, for readability. To get
an overview over all possible parameters look at
`polychromosims/params_default.py`
