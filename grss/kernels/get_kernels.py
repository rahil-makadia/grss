"""Download the custom and generic SPICE kernels from the
GRSS GitHub repository and the NAIF FTP server"""
import os

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

NAIF_SITE = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels'
SSD_SITE = 'https://ssd.jpl.nasa.gov/ftp'

# get the spice kernels if they are not already present
# de430 planets + de431 big16
if (not os.path.exists(f'{script_dir}/de430.bsp') or
    not os.path.exists(f'{script_dir}/sb431-n16s.bsp')):
    print('Downloading DE430/431 planet and big16 asteroid kernels...')
    os.system((f'curl --silent --show-error -o {script_dir}/de430.bsp '
                f'{NAIF_SITE}/spk/planets/de430.bsp'))
    os.system((f'curl --silent --show-error -o {script_dir}/sb431-n16s.bsp '
                f'{SSD_SITE}/xfr/sb431-n16s.bsp'))
# de440 planets + de441 big16
if (not os.path.exists(f'{script_dir}/de440.bsp') or
    not os.path.exists(f'{script_dir}/sb441-n16s.bsp')):
    print('Downloading DE440/441 planet and big16 asteroid kernels...')
    os.system((f'curl --silent --show-error -o {script_dir}/de440.bsp '
                f'{NAIF_SITE}/spk/planets/de440.bsp'))
    os.system((f'curl --silent --show-error -o {script_dir}/sb441-n16s.bsp '
                f'{SSD_SITE}/xfr/sb441-n16s.bsp'))

# get the earth orientation binary spice kernels and their comments if they are not already present
# latest earth pck
if not os.path.exists(f'{script_dir}/earth_latest_high_prec.bpc'):
    print('Downloading latest Earth binary PCK...')
    os.system((f'curl --silent --show-error -o {script_dir}/earth_latest_high_prec.bpc '
                f'{NAIF_SITE}/pck/earth_latest_high_prec.bpc'))
# historical earth pck
if not os.path.exists(f'{script_dir}/earth_720101_230601.bpc'):
    print('Downloading historical Earth binary PCK...')
    os.system((f'curl --silent --show-error -o {script_dir}/earth_720101_230601.bpc '
                f'{NAIF_SITE}/pck/earth_720101_230601.bpc'))
# predicted earth pck
if not os.path.exists(f'{script_dir}/earth_200101_990825_predict.bpc'):
    print('Downloading predicted Earth binary PCK...')
    os.system((f'curl --silent --show-error -o {script_dir}/earth_200101_990825_predict.bpc '
                f'{NAIF_SITE}/pck/earth_200101_990825_predict.bpc'))
# generic frame kernels
if not os.path.exists(f'{script_dir}/pck00011.tpc'):
    print('Downloading generic frame kernels...')
    os.system((f'curl --silent --show-error -o {script_dir}/pck00011.tpc '
                f'{NAIF_SITE}/pck/pck00011.tpc'))
