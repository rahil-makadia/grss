"""Download the debiasing data files from the JPL Solar System Dynamics Group's FTP server"""
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import requests

def check_and_download(url, local_path, note):
    """Download a bias-data archive to an existing directory, replacing its file.

    ``url`` is checked with HEAD before GET; HTTP/network failures raise rather
    than writing an error response into the destination. Returns None.
    """
    try:
        response = requests.head(url, timeout=30, allow_redirects=True)
        if response.status_code == 200:
            print(f'Downloading {note}...')
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            with open(local_path, 'wb') as f:
                f.write(response.content)
        else:
            msg = f'File {url} ({note}) does not exist on server. Code: {response.status_code}'
            raise FileNotFoundError(msg)
    except requests.RequestException as e:
        raise ConnectionError(f'Error checking/downloading {note}') from e

cwd = os.path.dirname(os.path.abspath(__file__))

FTP_SITE = 'https://ssd.jpl.nasa.gov/ftp/ssd/debias'
__all__ = ['check_and_download', 'download_debiasing_data', 'extract_debiasing_archive']


def extract_debiasing_archive(archive_path, destination):
    """Unpack regular files/directories, rejecting traversal and archive links.

    Extraction is staged before replacing same-named destination files. Other
    destination files are preserved. The input archive is not removed.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as staging:
        root = Path(staging).resolve()
        with tarfile.open(archive_path, 'r:gz') as archive:
            members = archive.getmembers()
            for member in members:
                target = (root/member.name).resolve()
                if (not target.is_relative_to(root)
                        or not (member.isfile() or member.isdir())):
                    raise ValueError(f'Unsafe debiasing archive member: {member.name}')
            for member in members:
                target = root/member.name
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.extractfile(member) as source, target.open('wb') as output:
                        shutil.copyfileobj(source, output)
        for source in root.rglob('*'):
            target = destination/source.relative_to(root)
            if not target.resolve().is_relative_to(destination.resolve()):
                raise ValueError(f'Unsafe debiasing destination: {target}')
            if source.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                os.replace(source, target)


def download_debiasing_data(destination=None):
    """Download and unpack missing low/high-resolution astrometric bias tables.

    Parameters
    ----------
    destination : str or os.PathLike, optional
        Output directory; defaults to the package's debias directory. Explicit
        calls perform network requests and replace files when unpacking an
        existing archive. Importing this module does not download or unpack.
    """
    cwd = os.fspath(destination) if destination is not None else os.path.dirname(__file__)
    os.makedirs(cwd, exist_ok=True)
    # get debiasing data files from the JPL Solar System Dynamics Group's FTP server
    # if lowres_data.tgz or lowres_data/ already exist, do not download
    if not os.path.exists(f'{cwd}/lowres_data.tgz') and not os.path.exists(f'{cwd}/lowres_data/'):
        check_and_download(f'{FTP_SITE}/debias_2018.tgz',
                            f'{cwd}/lowres_data.tgz', 'low-resolution debiasing data')
    if not os.path.exists(f'{cwd}/hires_data.tgz') and not os.path.exists(f'{cwd}/hires_data/'):
        check_and_download(f'{FTP_SITE}/debias_hires2018.tgz',
                            f'{cwd}/hires_data.tgz', 'high-resolution debiasing data')

    # Unpack only explicitly supplied/downloaded archives; preserve unrelated files.
    if os.path.exists(f'{cwd}/lowres_data.tgz'):
        extract_debiasing_archive(f'{cwd}/lowres_data.tgz', f'{cwd}/lowres_data')
        os.remove(f'{cwd}/lowres_data.tgz')
    if os.path.exists(f'{cwd}/hires_data.tgz'):
        extract_debiasing_archive(f'{cwd}/hires_data.tgz', f'{cwd}/hires_data')
        os.remove(f'{cwd}/hires_data.tgz')


if __name__ == '__main__':
    download_debiasing_data()
