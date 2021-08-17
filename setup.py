# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

"""PyGAMer: Geometry-preserving Adaptive Mesher

PyGAMer is a wrapper around the `GAMer <https://github.com/ctlee/gamer/`__ C++ library.
It enables users to perform mesh generation and manipulation tasks via Python scripts.
"""

import os
import sys
import subprocess
import re

def git_version():
    """Get the version from git describe

    Returns:
        string: version string or None if invalid
    """
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ["SYSTEMROOT", "PATH", "HOME"]:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env["LANGUAGE"] = "C"
        env["LANG"] = "C"
        env["LC_ALL"] = "C"
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, env=env)
        return out

    try:
        out = _minimal_ext_cmd(["git", "describe", "--tags", "--dirty", "--always"])
        GIT_REVISION = out.strip().decode("ascii")
    except (subprocess.SubprocessError, OSError):
        GIT_REVISION = None

    return GIT_REVISION

def standardize_version(version_string):
    """Standardize the version string
    
    Args:
        version_string (str): input version string 
    Returns:
        str: standardized version
    """
    VERSION_PATTERN = r"""
        v?
        (?:
            (?:(?P<epoch>[0-9]+)!)?                           # epoch
            (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
            (?P<pre>                                          # pre-release
                [-_\.]?
                (?P<pre_l>(alpha|beta|a|b|c|rc))
                (?P<pre_n>[0-9]+)?
            )?
            (?P<dev>                                          # dev release
                [-_\.]?
                (?P<dev_l>dev)
                [-_\.]?
                (?P<dev_n>[0-9]+)?
            )?
            (?P<meta>
                [-_\.]?
                (?P<commits_since>[0-9]+)?
                [-_\.]?
                (?P<sha>[a-z0-9]*)?
                [-_\.]? 
                (?P<dirty>dirty)?
            )?
        )
    """

    _regex = re.compile(
        r"^\s*" + VERSION_PATTERN + r"\s*$",
        re.VERBOSE | re.IGNORECASE,
    )

    match = _regex.match(version_string)
    if match:
        if match.group("release"):
            version = match.group("release")
            if match.group("pre_l"):
                version += match.group("pre_l")
                if match.group("pre_n"):
                    version += match.group("pre_n")
                else:
                    version += "0"
            if match.group("dev"):
                version += "dev"
                if match.group("dev_n"):
                    version += match.group("dev_n")
                else:
                    version += "0"
    else:
        version = "0.0.0"
    return version

version = git_version()
if version is None:
    # Git describe failed... read version from file
    with open("VERSION", "r") as f:
        version = f.readline()
    version = standardize_version(version)
else:
    version = standardize_version(version)

cmake_args = ["-DBUILD_PYGAMER=ON"]

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Natural Language :: English
Operating System :: OS Independent
Programming Language :: C++
Programming Language :: Python :: 3 :: Only
Programming Language :: Python :: Implementation :: CPython
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Mathematics
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Visualization
"""

# If building on readthedocs.io build with unix makefiles
_on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if _on_rtd:
    sys.argv.extend(["-G", "Unix Makefiles"])
    cmake_args.append("-DGAMER_DOCS=ON")

try:
    from skbuild import setup
except ImportError:
    print("\nERROR: scikit-build is required to build from source.", file=sys.stderr)
    print("Please run:", file=sys.stderr)
    print("", file=sys.stderr)
    print("  python -m pip install scikit-build\n", file=sys.stderr)
    print("  -- or --\n", file=sys.stderr)
    print("  conda install scikit-build", file=sys.stderr)
    sys.exit(1)

docs_require = [
    "sphinx",
    "breathe",
    "exhale",
    "sphinx-issues",
    "nbsphinx",
    "jupyter_sphinx",
    "threevis",
]
tests_require = ["pytest"]

setup(
    name="pygamer",
    version=version,
    maintainer="Christopher T. Lee",
    maintainer_email="ctlee@ucsd.edu",
    author="The GAMer Team",
    author_email="ctlee@ucsd.edu",
    url="https://github.com/ctlee/gamer",
    license="LGPLv2+",
    packages=["pygamer"],
    description=DOCLINES[0],
    long_description=open("README.md", encoding="utf8").read(),
    long_description_content_type="text/markdown",
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    classifiers=[c for c in CLASSIFIERS.split("\n") if c],
    keywords="meshing ",
    cmake_args=cmake_args,
    setup_requires=[
        "setuptools",
        "wheel",
        "scikit-build >= 0.10.0",
        "pytest-runner",
        "cmake >= 3.11",
    ],
    install_requires=["numpy>=1.8.0"],
    extras_require={
        "docs": docs_require,
        "test": tests_require,
    },
    zip_safe=False,
)
