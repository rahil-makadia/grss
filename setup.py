""" Setup script for grss package. """
# taken and adapted from setup.py at https://github.com/pybind/cmake_example
import os
import subprocess
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
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
        os.system(f"cp {ext.sourcedir}/build/prop_simulation* {self.build_lib}/grss/prop/")
        return

# get version from version.txt
with open("grss/version.txt", "r", encoding="utf-8") as f:
    ver = f.read().strip()
# run get_cspice in the extern folder if installing from source
if (not os.path.exists("./extern/cspice/lib/cspice.a") or
    not os.path.exists("./extern/cspice/include/SpiceUsr.h") ):
    subprocess.run(["python", "get_cspice.py"], cwd="./extern", check=True)
setup(
    version=ver,
    packages=["grss", "grss.fit", "grss.prop"],
    ext_modules=[CMakeExtension("prop_simulation")],
    cmdclass={"build_ext": CMakeBuild},
)
