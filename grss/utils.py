"""GRSS Python utilities subpackage"""
import os
import time
import gzip
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

def _download_ades_files(file):
    """
    Download a text file of the observatory code info from the Minor Planet Center.

    Parameters
    ----------
    file : str
        The name of the file to download from the Minor Planet Center.
    """
    print(f"Downloading ADES file {file} from the MPC...")
    files = {
        # "modes.json": "https://www.minorplanetcenter.net/iau/info/mode.json",
        # "catalogs.json": "https://www.minorplanetcenter.net/iau/info/astCat_photCat.json",
        "codes.json": "https://www.minorplanetcenter.net/Extended_Files/obscodes_extended.json.gz",
    }
    outname = f'{grss_path}/fit/{file}'
    if file == "codes.json":
        outname += ".gz"
    cmd = f'curl --silent --show-error -o {outname} {files[file]}'
    os.system(cmd)
    if file == "codes.json":
        with gzip.open(f'{outname}', 'rb') as f_in:
            with open(f'{grss_path}/fit/codes.json', 'w', encoding='utf-8') as f_out:
                f_out.write(f_in.read().decode('utf-8'))
        os.remove(f'{outname}')
    return None

def check_and_get_ades_files():
    """
    Download ADES JSON files from the Minor Planet Center.
    """
    # check if all files exist
    # files = ["modes.json", "catalogs.json", "codes.json"]
    files = ["codes.json"]
    for file in files:
        if not os.path.exists(f'{grss_path}/fit/{file}'):
            _download_ades_files(file)
        else:
            # check if file is older than 1 day
            file_age = time.time() - os.path.getmtime(f'{grss_path}/fit/{file}')
            if file_age > 86400:
                _download_ades_files(file)
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
    check_and_get_ades_files()
    # set openmp environment variable to number of cores
    os.environ['OMP_NUM_THREADS'] = str(os.cpu_count())
    return None
