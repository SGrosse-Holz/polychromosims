#!/home/simongh/anaconda3/bin/python3.7
# This script scans all the files in the main folder and checks that all
# occurences of params.(something) have a default value in params_defaults.py
# NOTE: this depends on that module always being called "params".

import sys,os
import re

from polychromosims import params_default

if __name__ == "__main__":
    has_default = [att for att in dir(params_default)]
    # Explicit exceptions:
    has_default += ["folder",
                    "py", # We occasionally talk about params.py in comments and strings
                    "gen_start_conf",
                    "cl_bondWiggleDist",
                    "cl_bondLength",
                    "extrusion_stepsPerBlock",
                    "extrusion_CTCFs",
                    "extrusion_totalRestarts",
                    "extrusion_blocksPerRestart"]

    # This should be the official spec for python variable names
    attrib_pattern = re.compile("params\.([a-zA-Z_][a-zA-Z0-9_]*)")
    
    # cd ../polychromosims :
    psims_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "polychromosims")

    missing = {} # Structure: {'attr_name' : [list of positions, a la 'file.py : line x']}
    for root, dirnames, filenames in os.walk(psims_dir):
        for filename in filenames:
            if not filename.endswith(".py"):
                continue
    
            with open(os.path.join(root, filename)) as curfile:
                lineno = 0
                for line in curfile:
                    lineno += 1
                    for match in attrib_pattern.finditer(line):
                        if not match[1] in has_default:
                            if not match[1] in missing.keys():
                                missing[match[1]] = []
                            missing[match[1]].append(os.path.relpath(os.path.join(root, filename),
                                                     psims_dir)+' : line {}'.format(lineno))

    if missing:
        print("found attributes with missing default values:")
        print("")
        for attr in missing.keys():
            print(attr, "in")
            for loc in missing[attr]:
                print(loc)
            print("")
    else:
        print("no missing defaults!")
