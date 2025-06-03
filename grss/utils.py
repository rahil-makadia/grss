"""GRSS Python utilities subpackage"""
import os
import time
import json
import requests

__all__ = [ 'default_kernel_path',
            'grss_path',
            'initialize'
]

grss_path = os.path.dirname(os.path.abspath(__file__))
default_kernel_path = f'{grss_path}/kernels/'

def _connected_to_internet(url='http://www.google.com/', timeout=30):
    try:
        _ = requests.head(url, timeout=timeout)
        return True
    except (requests.ConnectionError, requests.Timeout):
        print("No internet connection available.")
    return False

def _download_codes_file():
    """
    Download the observatory codes file from the Minor Planet Center.
    """
    print("Updating observatory codes from the MPC...")
    url = "https://data.minorplanetcenter.net/api/obscodes"
    response = requests.get(url, json={}, timeout=30)
    if not response.ok:
        raise ValueError("_download_codes_file Error: ", response.status_code, response.content)
    all_data = response.json()
    mpc_info_dict = {
        key: (
            float(val['longitude']),
            float(val['rhocosphi']),
            float(val['rhosinphi']),
        )
        for key, val in all_data.items()
            if all_data[key]['observations_type'] not in {'satellite', 'roving'} and key != '275'
    }
    fpath = f'{grss_path}/fit/codes.json'
    # write mpc_info_dict to file
    with open(fpath, 'w', encoding='utf-8') as f:
        json.dump(mpc_info_dict, f, indent=4)
    return None

def check_and_get_codes_file():
    """
    Check if the observatory codes file exists and is up to date. If not, download it.
    """
    fpath = f'{grss_path}/fit/codes.json'
    if not os.path.exists(fpath):
        _download_codes_file()
    else:
        # check if file is older than 1 day
        file_age = time.time() - os.path.getmtime(fpath)
        if file_age > 86400:
            _download_codes_file()
    return None

def initialize():
    """
    GRSS library initialization function. This function will download the
    necessary data files for the GRSS library to work that aren't included
    in the package distribution (e.g. SPICE kernels, debiasing data).

    Returns
    -------
    None : NoneType
        None
    """
    internet = _connected_to_internet()
    if not internet:
        print("Please connect to the internet to download/update the necessary data files.")
        return None
    # run the get_debiasing_data.py script
    os.system(f'python {grss_path}/debias/get_debiasing_data.py')
    # run the get_kernels.py script
    os.system(f'python {grss_path}/kernels/get_kernels.py')
    # get the observatory codes file
    check_and_get_codes_file()
    # set openmp environment variable to number of cores
    os.environ['OMP_NUM_THREADS'] = str(os.cpu_count())
    return None
