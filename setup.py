""" Setup script for grss package. """
# taken and adapted from setup.py at https://github.com/pybind/cmake_example
import os
import subprocess
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    """
    CMake extension for setuptools.

    Parameters
    ----------
    Extension : setuptools.Extension
        Extension object from setuptools.
    """
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())

class CMakeBuild(build_ext):
    """
    CMake build class for setuptools.

    Parameters
    ----------
    build_ext : setuptools.command.build_ext.build_ext
        Build extension object from setuptools.
    """

    def build_extension(self, ext: CMakeExtension) -> None:
        subprocess.run(["./build_cpp.sh"], cwd=ext.sourcedir, check=True)
        binary_created = [f for f in os.listdir(f"{ext.sourcedir}/build/")
                            if f.startswith("libgrss")]
        if not binary_created:
            raise FileNotFoundError("libgrss binary for C++ source code not found "
                                    "in cmake build directory")
        os.system(f"cp {ext.sourcedir}/build/libgrss* {self.build_lib}/grss/")
        return

# run get_cspice in the extern folder if installing from source
if (not os.path.exists("./extern/cspice/lib/cspice.a") or
    not os.path.exists("./extern/cspice/include/SpiceUsr.h") ):
    subprocess.run(["python", "get_cspice.py"], cwd="./extern", check=True)

setup(
    ext_modules=[CMakeExtension("libgrss")],
    cmdclass={"build_ext": CMakeBuild}
)
