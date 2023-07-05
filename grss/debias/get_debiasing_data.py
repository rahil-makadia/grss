"""Download the debiasing data files from the JPL Solar System Dynamics Group's FTP server"""
import os

cwd = os.path.dirname(os.path.abspath(__file__))

FTP_SITE = 'ftp://ssd.jpl.nasa.gov/pub/ssd/debias'
# get debiasing data files from the JPL Solar System Dynamics Group's FTP server
# if lowres_data.tgz or lowres_data/ already exist, do not download
if not os.path.exists(f'{cwd}/lowres_data.tgz') and not os.path.exists(f'{cwd}/lowres_data/'):
    FTP_FILE = 'debias_2018.tgz'
    os.system(f'wget --no-verbose --no-clobber {FTP_SITE}/{FTP_FILE} -O {cwd}/lowres_data.tgz')
if not os.path.exists(f'{cwd}/hires_data.tgz') and not os.path.exists(f'{cwd}/hires_data/'):
    FTP_FILE = 'debias_hires2018.tgz'
    os.system(f'wget --no-verbose --no-clobber {FTP_SITE}/{FTP_FILE} -O {cwd}/hires_data.tgz')

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
