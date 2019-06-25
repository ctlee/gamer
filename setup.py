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

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH', 'HOME']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, env=env)
        return out

    try:
        out = _minimal_ext_cmd(['git', 'describe', '--tags', '--dirty=.dirty', '--always'])
        GIT_REVISION = out.strip().decode('ascii')
    except (subprocess.SubprocessError, OSError):
        GIT_REVISION = "Unknown"

    return GIT_REVISION

version = git_version()
if not version == "Unknown":
    match = re.search('v(\d+)\.(\d+)\.(\d+)-*(alpha|beta|dev|)-*([A-Za-z0-9_-]*)\.*(dirty|)', version)
    if match:
        version = "%s.%s.%s"%(match.group(1), match.group(2), match.group(3))
    else:
        version = "0.0.0"
else:
    with open('VERSION', 'r') as f:
       version = f.readline()
    match = re.search('(\d+)\.(\d+)\.(\d+)-*(alpha|beta|dev|)', version)
    if match:
        version = "%s.%s.%s"%(match.group(1), match.group(2), match.group(3))
    else:
        version = "0.0.0"

cmake_args=['-DBUILD_PYGAMER=ON']

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
_on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if(_on_rtd):
    sys.argv.extend(['-G','Unix Makefiles'])
    cmake_args.append('-DGAMER_DOCS=ON')

from skbuild import setup

docs_require = ["sphinx", "breathe", "exhale", "sphinx-issues", "nbsphinx",
"jupyter_sphinx", "threevis"]
tests_require = ["pytest"]

setup(
    name='pygamer',
    version=version,
    maintainer='Christopher T. Lee',
    maintainer_email='ctlee@ucsd.edu',
    author='The GAMer Team',
    author_email='ctlee@ucsd.edu',
    url='https://github.com/ctlee/gamer',
    license='LGPLv2+',
    packages=["pygamer"],
    description=DOCLINES[0],
    long_description=open('README.md', encoding='utf8').read(),
    long_description_content_type="text/markdown",
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    classifiers=[c for c in CLASSIFIERS.split('\n') if c],
    keywords='meshing ',
    cmake_args=cmake_args,
    setup_requires=["setuptools", "wheel", "scikit-build >= 0.10.0", "pytest-runner", "cmake >= 3.11"],
    install_requires=["numpy>=1.8.0"],
    extras_require={
        "docs" : docs_require,
        "test" : tests_require,
    },
    zip_safe=False,
)