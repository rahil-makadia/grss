# from the spiceypy cspice downloader
"""Download the CSPICE package from the NAIF FTP server and unpack it."""
import platform
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

CSPICE_SRC_DIR = "CSPICE_SRC_DIR"
CSPICE_SHARED_LIB = "CSPICE_SHARED_LIB"
CSPICE_NO_PATCH = "CSPICE_NO_PATCH"

host_OS = platform.system()
host_arch = platform.machine()
# Check if platform is supported
os_supported = host_OS in ("Linux", "Darwin")
# Get platform is Unix-like OS or not
is_unix = host_OS in ("Linux", "Darwin")
# Get current working directory
ROOT_DIR = str(Path(os.path.realpath(__file__)).parent)
# Make the directory path for cspice
CSPICE_DIR = os.environ.get(CSPICE_SRC_DIR, os.path.join(ROOT_DIR, "cspice"))
# if we need to cross compile or compile for arm64
is_macos_arm = host_OS == "Darwin" and (
    host_arch == "arm64" or os.environ.get("ARCHFLAGS", "") == "-arch arm64"
)
# versions
SPICE_VERSION = "N0067"
SPICE_NUM_V = "67"


class GetCSPICE(object):
    """
    Class that support the download from the NAIF FTP server of the required
    CSPICE package for the architecture used by the Python distribution that
    invokes this module.  By default the CSPICE Toolkit version N0066 is
    downloaded and unpacked on the directory where this module is located.

    Arguments
    ---------
    :argument version: String indicating the required version of the CSPICE
                       Toolkit. By default it is 'N0066'.
    :type: str

    """

    # This class variable will be used to store the CSPICE package in memory.
    _local = None

    # Supported distributions
    _dists = {
        # system   arch        distribution name           extension
        # -------- ----------  -------------------------   ---------
        ("Darwin", "x86_64", "64bit"): ("MacIntel_OSX_AppleC_64bit", "tar.Z"),
        ("Darwin", "arm64", "64bit"): ("MacM1_OSX_clang_64bit", "tar.Z"),
        ("Linux", "x86_64", "64bit"): ("PC_Linux_GCC_64bit", "tar.Z"),
        ("Linux", "aarch64", "64bit"): ("PC_Linux_GCC_64bit", "tar.Z"),
    }

    def __init__(self, version=SPICE_VERSION, dst=None):
        """Init method that uses either the default N0066 toolkit version token
        or a user provided one.
        """
        try:
            # Get the remote file path for the Python architecture that
            # executes the script.
            distribution, self._ext = self._distribution_info()
        except KeyError:
            print("GRSS currently does not support your system.")
        else:
            cspice = f"cspice.{self._ext}"
            self._rcspice = (
                f"https://naif.jpl.nasa.gov/pub/naif/misc"
                f"/toolkit_{version}/C/{distribution}/packages"
                f"/{cspice}"
            )

            # Setup the local directory (where the package will be downloaded)
            if dst is None:
                dst = os.path.realpath(os.path.dirname(__file__))
            self._root = dst

            # Download the file
            print(f"Downloading CSPICE for {distribution}...")
            attempts = 10  # Let's try a maximum of attempts for getting SPICE
            while attempts:
                attempts -= 1
                try:
                    self._download()
                except RuntimeError as error:
                    print(
                        f"Download failed with URLError: {error}, trying again after "
                        "15 seconds!"
                    )
                    time.sleep(15)
                else:
                    # Unpack the file
                    print("Unpacking... (this may take some time!)")
                    self._unpack()
                    # We are done.  Let's return to the calling code.
                    break
        return

    def _distribution_info(self):
        """Creates the distribution name and the expected extension for the
        CSPICE package and returns it.

        :return (distribution, extension) tuple where distribution is the best
                guess from the strings available within the platform_urls list
                of strings, and extension is either "zip" or "tar.Z" depending
                on whether we are dealing with a Windows platform or else.
        :rtype: tuple (str, str)

        :raises: KeyError if the (system, machine) tuple does not correspond
                    to any of the supported GRSS environments.
        """
        print("Gathering information...")
        system = platform.system()
        processor = platform.processor()
        machine = platform.machine()
        # Cygwin system is CYGWIN-NT-xxx.
        system = "cygwin" if "CYGWIN" in system else system
        cpu_bits = "64bit" if sys.maxsize > 2 ** 32 else "32bit"
        if machine in ("x86", "x86_64", "AMD64", "i686"):
            machine = "x86_64"
        if is_macos_arm:
            print("either running on apple arm64 or cross-compiling")
            machine = "arm64"
        print("SYSTEM:   ", system)
        print("PROCESSOR:", processor)
        print("MACHINE:  ", machine, cpu_bits)
        if machine in ("i386", "x86_32") or cpu_bits == "32bit":
            raise RuntimeError("32bit bit builds are not supported")
        return self._dists[(system, machine, cpu_bits)]

    def _download(self):
        """Support function that encapsulates the OpenSSL transfer of the CSPICE
        package to the self._local io.ByteIO stream.

        :raises RuntimeError if there has been any issue with the HTTPS
                communication

        .. note::

            Handling of CSPICE downloads from HTTPS
            ---------------------------------------
            Some Python distributions may be linked to an old version of OpenSSL
            which will not let you connect to NAIF server due to recent SSL cert
            upgrades on the JPL servers.  Moreover, versions older than
            OpenSSL 1.0.1g are known to contain the 'the Heartbleed Bug'.
            Therefore this method provides two different implementations for the
            HTTPS GET call to the NAIF server to download the required CSPICE
            distribution package.
        """
        cmd = f"wget --no-clobber {self._rcspice}"
        subprocess.run(cmd,shell=True,check=True)
        return

    def _unpack(self):
        """Unpacks the CSPICE package on the given root directory. Note that
        Package could either be the zipfile.ZipFile class for Windows platforms
        or tarfile.TarFile for other platforms.
        """
        cspice = f"cspice.{self._ext}"
        opts = '-xzf'
        cmd = f"tar {opts} {cspice} -C {self._root}"
        subprocess.run(cmd,shell=True,check=True)
        os.remove(cspice)
        return


def copy_supplements() -> None:
    """
    Copy supplement files (patches, windows build files)
    to cspice directory
    """
    cwd = os.getcwd()
    patches = list(Path().cwd().glob("*.patch"))
    print("copy supplements to: ", CSPICE_DIR, flush=True)
    for patch in patches:
        shutil.copy(patch, CSPICE_DIR)
    if host_OS == "Windows":
        windows_files = [Path("./makeDynamicSpice.bat"), Path("./cspice.def")]
        cspice_src_dir = os.path.join(CSPICE_DIR, "cspice", "src", "cspice")
        for win_file in windows_files:
            shutil.copy(win_file, cspice_src_dir)
    os.chdir(cwd)
    return


def apply_patches() -> None:
    """
    Apply patches to cspice source code
    """
    if int(SPICE_NUM_V) != 66:
        return
    cwd = os.getcwd()
    os.chdir(CSPICE_DIR)
    iswin = "-windows" if host_OS == "Windows" else ""
    patches = [
        f"0001-patch-for-n66-dskx02.c{iswin}.patch",
        f"0002-patch-for-n66-subpnt.c{iswin}.patch",
    ]
    if is_macos_arm:
        patches.append("0004_inquire_unistd.patch")
    for patch in patches:
        try:
            print(f"Applying Patch {patch}", flush=True)
            subprocess.run(["git", "apply", "--reject", patch], check=True)
        except subprocess.CalledProcessError as cpe:
            raise cpe
    os.chdir(cwd)
    return


def prepare_cspice() -> None:
    """
    Prepare temporary cspice source directory,
    If not provided by the user, or if not readable, download a fresh copy
    :return: None
    """
    cwd = os.getcwd()
    cspice_root_dir = str((Path(ROOT_DIR)).absolute())
    tmp_cspice_src_dir = os.path.join(cspice_root_dir, "cspice")
    if os.access(CSPICE_DIR, os.R_OK):
        # remove any existing cspice dir
        shutil.rmtree(CSPICE_DIR)
    else:
        print("Downloading CSPICE src from NAIF", flush=True)
    Path(tmp_cspice_src_dir).mkdir(exist_ok=True, parents=True)
    GetCSPICE(dst=cspice_root_dir)
    os.chdir(cwd)
    # okay now copy any and all files needed for building
    return


def build_cspice():
    """
    Builds cspice
    :return: absolute path to new compiled shared library and obsolute path to header files
    """
    cwd = os.getcwd()
    if is_unix:
        minv_flag = ""
        libname = "libcspice.a"
        target = "-target arm64-apple-macos11" if is_macos_arm else ""
        if host_OS == "Darwin":
            # extra_flags = f"-install_name @rpath/{libname}"
            extra_flags = "-dynamiclib"
            minv_flag = "-mmacosx-version-min=10.15"
        else:
            extra_flags = f"-Wl,{libname}"
        destination = CSPICE_DIR
        os.chdir(destination)
        cmds = [
            f"gcc {target} {minv_flag} -Iinclude -c -fPIC -O2 -ansi -w ./cspice/src/cspice/*.c",
            f"gcc {target} {extra_flags} {minv_flag} -fPIC -O2 -lm *.o -o {libname}",
        ]
    elif host_OS == "Windows":
        destination = os.path.join(CSPICE_DIR, "cspice", "src", "cspice")
        os.chdir(destination)
        cmds = ["makeDynamicSpice.bat"]
    else:
        os.chdir(cwd)
        raise NotImplementedError(f"non implemented host os for build {host_OS}")
    try:
        for cmd in cmds:
            _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        os.chdir(cwd)
    # get the built  library
    shared_lib_path = [
        str(p.absolute())
        for p in Path(destination).glob("*.*")
        if p.suffix in (".lib", ".a")
    ]
    if len(shared_lib_path) != 1:
        os.chdir(cwd)
        raise RuntimeError(('Could not find built static library'
                                f'in {list(Path(destination).glob("*.*"))}'))
    shared_lib_path = shared_lib_path[0]
    print(shared_lib_path, flush=True)
    #Get the include path from the cspice src
    shared_include_path = os.path.join(destination,"cspice","include")
    print(shared_include_path, flush=True)
    os.chdir(cwd)
    return [shared_lib_path,shared_include_path]


def main() -> None:
    """
    Main routine to build or not build cspice
    expected tmp src dir layout
    /tmpdir/ < top level temporary directory from TemporaryDirectory
        /cspice/ < directory containing contents of source code for cspice
            /lib
            /bin
            /src/cspice/
            ...
    :return: None
    """
    cwd = os.getcwd()
    # set final destination for cspice dynamic library
    destination = os.path.join(
        ROOT_DIR,
        "extern",
        "cspice",
        "lib",
        "cspice.a" if is_unix else "cspice.lib",
    )
    # check if the shared library already exists, if it does we are done
    if Path(destination).is_file():
        print(
            "Done! shared library for cspice already exists in destination. Done!",
            flush=True,
        )
        return
    # next see if cspice shared library is provided
    shared_library_path = os.environ.get(CSPICE_SHARED_LIB)
    if shared_library_path is not None:
        print("User has provided a shared library...", flush=True)
        # what if we can't read the file? we need to jump to building it... doubt this happens
    else:
        # okay, now we either are given a src dir, have already downloaded it, or don't have it
        print("Preparing cspice", flush=True)
        prepare_cspice()
        # # add the patches
        # print("Copying supplements", flush=True)
        # copy_supplements()
        # # okay now that we have the source in a writeable directory we can apply patches
        # if CSPICE_NO_PATCH not in os.environ:
        #     print("Apply patches", flush=True)
        #     apply_patches()
        # now build
        # print("Building cspice", flush=True)
        # shared_library_path,shared_include_path = build_cspice()

    os.chdir(cwd)
    print("Done!")

if __name__ == "__main__":
    main()
