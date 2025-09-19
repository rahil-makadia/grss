"""Download the debiasing data files from the JPL Solar System Dynamics Group's FTP server"""
import os
import requests

def check_and_download(url, local_path, note):
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

cwd = os.path.dirname(os.path.abspath(__file__))

FTP_SITE = 'https://ssd.jpl.nasa.gov/ftp/ssd/debias'
# get debiasing data files from the JPL Solar System Dynamics Group's FTP server
# if lowres_data.tgz or lowres_data/ already exist, do not download
if not os.path.exists(f'{cwd}/lowres_data.tgz') and not os.path.exists(f'{cwd}/lowres_data/'):
    check_and_download(f'{FTP_SITE}/debias_2018.tgz',
                        f'{cwd}/lowres_data.tgz', 'low-resolution debiasing data')
if not os.path.exists(f'{cwd}/hires_data.tgz') and not os.path.exists(f'{cwd}/hires_data/'):
    check_and_download(f'{FTP_SITE}/debias_hires2018.tgz',
                        f'{cwd}/hires_data.tgz', 'high-resolution debiasing data')

# extract the data files
# create the directories if they don't exist, delete the files if they do
# if tarballs don't exist, do nothing
if os.path.exists(f'{cwd}/lowres_data.tgz'):
    os.system(f'mkdir -p {cwd}/lowres_data')
    os.system(f'rm -f {cwd}/lowres_data/*')
    os.system(f'tar -xzf {cwd}/lowres_data.tgz -C {cwd}/lowres_data')
if os.path.exists(f'{cwd}/hires_data.tgz'):
    os.system(f'mkdir -p {cwd}/hires_data')
    os.system(f'rm -f {cwd}/hires_data/*')
    os.system(f'tar -xzf {cwd}/hires_data.tgz -C {cwd}/hires_data')

# remove the tarballs
os.system(f'rm -f {cwd}/lowres_data.tgz')
os.system(f'rm -f {cwd}/hires_data.tgz')
