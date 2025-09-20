"""Download the custom and generic SPICE kernels from the
GRSS GitHub repository and the NAIF FTP server"""
import os
import time
import requests

def check_and_download(url, local_path, note, force=False):
    if os.path.exists(local_path) and not force:
        return
    try:
        response = requests.head(url, timeout=30, allow_redirects=True)
        if response.status_code == 200:
            print(f'Downloading {note}...')
            response = requests.get(url, timeout=30)
            with open(local_path, 'wb') as f:
                f.write(response.content)
        else:
            msg = f'File {url} ({note}) does not exist on server. Code: {response.status_code}'
            raise FileNotFoundError(msg)
    except requests.RequestException as e:
        raise ConnectionError(f'Error checking/downloading {note}') from e

# get the path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))

NAIF_SITE = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels'
SSD_SITE = 'https://ssd.jpl.nasa.gov/ftp'

# get the spice kernels if they are not already present
# de430 planets + de431 big16
check_and_download(f'{NAIF_SITE}/spk/planets/de430.bsp',
                    f'{script_dir}/de430.bsp', 'DE430 planet kernel')
check_and_download(f'{SSD_SITE}/xfr/sb431-n16s.bsp',
                    f'{script_dir}/sb431-n16s.bsp', 'DE431 big16 asteroid kernel')
# de440/441 planets + de441 big16
check_and_download(f'{SSD_SITE}/eph/planets/bsp/de440.bsp',
                    f'{script_dir}/de440.bsp', 'DE440 planet kernel')
check_and_download(f'{SSD_SITE}/xfr/sb441-n16s.bsp',
                    f'{script_dir}/sb441-n16s.bsp', 'short-term DE441 big16 asteroid kernel')
check_and_download(f'{SSD_SITE}/eph/planets/bsp/de441.bsp',
                    f'{script_dir}/de441.bsp', 'DE441 planet kernel')
check_and_download(f'{SSD_SITE}/eph/small_bodies/asteroids_de441/sb441-n16.bsp',
                    f'{script_dir}/sb441-n16.bsp', 'long-term DE441 big16 asteroid kernel')

# latest earth pck
latest_earth_pck_exists = os.path.exists(f'{script_dir}/earth_latest.bpc')
latest_earth_pck_old = (latest_earth_pck_exists and
                        time.time() - os.path.getmtime(f'{script_dir}/earth_latest.bpc')
                        > 86400)
if not latest_earth_pck_exists or latest_earth_pck_old:
    check_and_download(f'{NAIF_SITE}/pck/earth_latest_high_prec.bpc',
                       f'{script_dir}/earth_latest.bpc', 'latest Earth binary PCK', force=True)
# historical earth pck
check_and_download(f'{NAIF_SITE}/pck/earth_620120_250826.bpc',
                    f'{script_dir}/earth_historic.bpc', 'historic Earth binary PCK')
# predicted earth pck
check_and_download(f'{NAIF_SITE}/pck/earth_2025_250826_2125_predict.bpc',
                    f'{script_dir}/earth_predict.bpc', 'predicted Earth binary PCK')
# moon pck
check_and_download(f'{NAIF_SITE}/pck/moon_pa_de440_200625.bpc',
                    f'{script_dir}/moon_pa_de440.bpc', 'Moon PCK for DE440')
# generic frame kernels
check_and_download(f'{NAIF_SITE}/pck/pck00011.tpc',
                    f'{script_dir}/pck00011.tpc', 'generic frame kernel')
