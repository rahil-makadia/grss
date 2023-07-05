"""Download the custom and generic SPICE kernels from the
GRSS GitHub repository and the NAIF FTP server"""
import os

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

GRSS_SITE = 'https://github.com/rahil-makadia/grss/raw/dev/grss/kernels'
NAIF_SITE = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels'
# get the custom spice kernels if they are not already present
# de431 planets + big16 1950-2350
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de431_1950_2350.bsp '
            f'-O {script_dir}/planets_big16_de431_1950_2350.bsp'))
# de441 planets + big16 1950-2350
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de441_1950_2350.bsp '
            f'-O {script_dir}/planets_big16_de441_1950_2350.bsp'))

# get the latest spice leap second kernel
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/lsk/latest_leapseconds.tls '
            f'-O {script_dir}/latest_leapseconds.tls'))
# get the earth orientation binary spice kernels and their comments if they are not already present
# latest earth pck
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_latest_high_prec.cmt '
            f'-O {script_dir}/earth_latest_high_prec.cmt'))
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_latest_high_prec.bpc '
            f'-O {script_dir}/earth_latest_high_prec.bpc'))
# historical earth pck
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_720101_230601.bpc '
            f'-O {script_dir}/earth_720101_230601.bpc'))
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_720101_230601.cmt '
            f'-O {script_dir}/earth_720101_230601.cmt'))
# predicted earth pck
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_200101_990825_predict.bpc '
            f'-O {script_dir}/earth_200101_990825_predict.bpc'))
os.system((f'wget --no-verbose --no-clobber {NAIF_SITE}/pck/earth_200101_990825_predict.cmt '
            f'-O {script_dir}/earth_200101_990825_predict.cmt'))

# open the spice meta-kernels and update the line that defines
# the PATH_VALUES variable to point to the same directory as this script
meta_kernels = [
    f'{script_dir}/planets_big16_de431_1950_2350.tm',
    f'{script_dir}/planets_big16_de441_1950_2350.tm'
]
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
