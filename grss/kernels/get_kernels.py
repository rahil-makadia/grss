"""Download the custom and generic SPICE kernels from the
GRSS GitHub repository and the NAIF FTP server"""
import os
import time

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

NAIF_SITE = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels'
SSD_SITE = 'https://ssd.jpl.nasa.gov/ftp'

# get the spice kernels if they are not already present
# de430 planets + de431 big16
mb430_exists = os.path.exists(f'{script_dir}/de430.bsp')
sb431_exists = os.path.exists(f'{script_dir}/sb431-n16s.bsp')
if (not mb430_exists or not sb431_exists):
    print('Downloading DE430/431 planet and big16 asteroid kernels...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/de430.bsp '
                f'{NAIF_SITE}/spk/planets/de430.bsp'))
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/sb431-n16s.bsp '
                f'{SSD_SITE}/xfr/sb431-n16s.bsp'))
# de440 planets + de441 big16
mb440_exists = os.path.exists(f'{script_dir}/de440.bsp')
sb441s_exists = os.path.exists(f'{script_dir}/sb441-n16s.bsp')
mb441_exists = os.path.exists(f'{script_dir}/de441.bsp')
sb441_exists = os.path.exists(f'{script_dir}/sb441-n16.bsp')
if (not mb440_exists or not mb441_exists or not sb441s_exists or not sb441_exists):
    print('Downloading DE440/441 planet and big16 asteroid kernels...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/de440.bsp '
                f'{SSD_SITE}/eph/planets/bsp/de440.bsp'))
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/sb441-n16s.bsp '
                f'{SSD_SITE}/xfr/sb441-n16s.bsp'))
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/de441.bsp '
                f'{SSD_SITE}/eph/planets/bsp/de441.bsp'))
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/sb441-n16.bsp '
                f'{SSD_SITE}/eph/small_bodies/asteroids_de441/sb441-n16.bsp'))

# get the earth orientation binary spice kernels and their comments if they are not already present
# latest earth pck
latest_earth_pck_exists = os.path.exists(f'{script_dir}/earth_latest.bpc')
latest_earth_pck_old = (latest_earth_pck_exists and
                        time.time() - os.path.getmtime(f'{script_dir}/earth_latest.bpc')
                        > 86400)
if not latest_earth_pck_exists or latest_earth_pck_old:
    print('Downloading latest Earth binary PCK...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/earth_latest.bpc '
                f'{NAIF_SITE}/pck/earth_latest_high_prec.bpc'))
# historical earth pck
historical_earth_pck_exists = os.path.exists(f'{script_dir}/earth_historic.bpc')
if not historical_earth_pck_exists:
    print('Downloading historic Earth binary PCK...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/earth_historic.bpc '
                f'{NAIF_SITE}/pck/earth_620120_240827.bpc'))
# predicted earth pck
predicted_earth_pck_exists = os.path.exists(f'{script_dir}/earth_predict.bpc')
if not predicted_earth_pck_exists:
    print('Downloading predicted Earth binary PCK...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/earth_predict.bpc '
                f'{NAIF_SITE}/pck/earth_200101_990827_predict.bpc'))
# moon pck
moon_pa_de440_exists = os.path.exists(f'{script_dir}/moon_pa_de440.bpc')
if not moon_pa_de440_exists:
    print('Downloading Moon PCK for DE440...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/moon_pa_de440.bpc '
                f'{NAIF_SITE}/pck/moon_pa_de440_200625.bpc'))
# generic frame kernels
generic_frame_kernel_exists = os.path.exists(f'{script_dir}/pck00011.tpc')
if not generic_frame_kernel_exists:
    print('Downloading generic frame kernel...')
    os.system((f'curl --connect-timeout 30 --show-error -o {script_dir}/pck00011.tpc '
                f'{NAIF_SITE}/pck/pck00011.tpc'))
