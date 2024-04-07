"""GRSS Python utilities subpackage"""
import os
import datetime
from requests import request

__all__ = [ 'default_kernel_path',
            'grss_path',
            'initialize'
]

grss_path = os.path.dirname(os.path.abspath(__file__))
default_kernel_path = f'{grss_path}/kernels/'

def _download_obs_codes_file():
    """
    Download a text file of the observatory code info from the Minor Planet Center.
    """
    print("Downloading observatory codes file from the MPC...")
    url = "https://www.minorplanetcenter.net/iau/lists/ObsCodes.html"
    req = request("GET", url, timeout=30)
    text = req.content.replace(b'<pre>',b'').replace(b'</pre>',b'')

    with open(f'{grss_path}/fit/obs_codes.txt', "wb") as file:
        now = datetime.datetime.now()
        file.write((f'{now.year}, {now.month}, {now.day}, '
                    f'{now.hour}, {now.minute}, {now.second}\n'.encode()))
        file.write(text[1:-2])
    return None

def get_obs_codes_file():
    # sourcery skip: extract-method
    """
    Download a text file of the observatory code info from the Minor Planet
    Center. If it already exists, update it if the last download was more than 1
    week ago.
    """
    # check if file exists
    if os.path.isfile(f'{grss_path}/fit/obs_codes.txt'):
        with open(f'{grss_path}/fit/obs_codes.txt', 'r', encoding='utf-8') as file:
            now = datetime.datetime.now()
            data = file.read()
            last_updated = data.split('\n')[0]
            last_updated = datetime.datetime(*[int(i) for i in last_updated.split(',')])
            diff_days = (now - last_updated).days
            if diff_days >= 1.0:
                _download_obs_codes_file()
    else:
        _download_obs_codes_file()
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
    # run the get_debiasing_data.py script
    os.system(f'python {grss_path}/debias/get_debiasing_data.py')
    # run the get_kernels.py script
    os.system(f'python {grss_path}/kernels/get_kernels.py')
    # get the observatory codes file
    get_obs_codes_file()
    # set openmp environment variable to number of cores
    os.environ['OMP_NUM_THREADS'] = str(os.cpu_count())
    return None
