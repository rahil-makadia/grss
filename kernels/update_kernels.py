# open the spice meta-kernels and update the line that defines 
# the PATH_VALUES variable to point to the same directory as this script

import os

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

# open the meta-kernel
meta_kernels = ['./planets_big16_de431_1950_2350.tm', './planets_big16_de441_1950_2350.tm']
for mk in meta_kernels:
    # read the meta-kernel and find the line that defines the PATH_VALUES variable
    with open(mk, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'PATH_VALUES' in line:
                # update the path to the directory containing this script
                lines[i] = f"    PATH_VALUES  = ( '{script_dir}" + "' )\n"
                break
    # write the updated meta-kernel
    with open(mk, 'w', encoding='utf-8') as f:
        f.writelines(lines)
