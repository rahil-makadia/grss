"""Download the custom and generic SPICE kernels from the
GRSS GitHub repository and the NAIF FTP server"""
import os
import sys

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

GRSS_SITE = 'https://github.com/rahil-makadia/grss/raw/dev/grss/kernels'
NAIF_SITE = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels'
# get the custom spice kernels if they are not already present
# de431 planets + big16 1950-2350
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de431_1950_2350.bsp '
            f'-O {script_dir}/planets_big16_de431_1950_2350.bsp'))
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de431_1950_2350.tm '
            f'-O {script_dir}/planets_big16_de431_1950_2350.tm'))
# de441 planets + big16 1950-2350
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de441_1950_2350.bsp '
            f'-O {script_dir}/planets_big16_de441_1950_2350.bsp'))
os.system((f'wget --no-verbose --no-clobber {GRSS_SITE}/planets_big16_de441_1950_2350.tm '
            f'-O {script_dir}/planets_big16_de441_1950_2350.tm'))

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

# run this code if the no-tm-overwrite flag argument is not present
if len(sys.argv) > 1:
    tm_overwrite = sys.argv[1] != '--no-tm-overwrite'
else:
    tm_overwrite = True
if tm_overwrite:
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
            num_chunks = 0
            for i, line in enumerate(lines):
                if 'PATH_VALUES' in line and 'placeholder' in line:
                    # update the path to the directory containing this script
                    # if script_dir is more that 40 characters long, then break it up
                    # into chunks of 40 characters each
                    cutoff = 40
                    if len(script_dir) > cutoff:
                        num_chunks, remainder = divmod(len(script_dir), cutoff)
                        chunks = [script_dir[i*cutoff:(i+1)*cutoff] for i in range(num_chunks)]
                        if remainder > 0:
                            chunks.append(script_dir[-remainder:])
                            num_chunks += 1
                        lines[i] = f"    PATH_VALUES  = ( '{chunks[0]}',\n"
                        for chunk in chunks[1:]:
                            end_char = " )" if chunk == chunks[-1] else ","
                            lines[i] += f"                     '{chunk}'{end_char}\n"
                    else:
                        num_chunks = 1
                        lines[i] = f"    PATH_VALUES  = ( '{script_dir}" + "' )\n"
                if 'PATH_SYMBOLS' in line and "'GRSS'" in line and num_chunks > 1:
                    # replace PATH_SYMBOLS = ( 'GRSS' ) with PATH_SYMBOLS = ( 'GRSS_1', 'GRSS_2', ... )
                    lines[i] = "    PATH_SYMBOLS = ( 'GRSS1',\n"
                    for j in range(2, num_chunks+1):
                        end_char = " )" if j == num_chunks else ","
                        lines[i] += f"                     'GRSS{j}'{end_char}\n"
                if '$GRSS' in line and num_chunks > 1:
                    # replace '$GRSS' with '$GRSS1$GRSS2$GRSS3...' according to the number of chunks
                    replacement_str = '$GRSS1'
                    for j in range(2, num_chunks+1):
                        replacement_str += f'$GRSS{j}'
                    lines[i] = line.replace('$GRSS', replacement_str)
        # write the updated meta-kernel
        with open(mk, 'w', encoding='utf-8') as f:
            f.writelines(lines)
