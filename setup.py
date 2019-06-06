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

PyGAMer is a wrapper around the `GAMer <https://github.com/ctlee/gamer/`_ C++ library.
It enables users to perform mesh generation and manipulation tasks via Python scripts.
"""

import os
import sys

cmake_args=['-DBUILD_PYGAMER=ON']

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv22+)
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

build_require = ["setuptools", "wheel", "scikit-build >= 0.10.0", "cmake >= 3.11", "ninja"]
docs_require = ["sphinx", "breathe", "exhale", "sphinx-issues"]
tests_require = ["pytest"]

setup(
    name='pygamer',
    version='2.0.2',
    maintainer='Christopher T. Lee',
    maintainer_email='ctlee@ucsd.edu',
    author='The GAMer Team',
    author_email='',
    url='https://github.com/cltee/gamer',
    license='LGPLv2+',
    description=DOCLINES[0],
    long_description='Python wrapper around the GAMer C++ library for mesh generation.',
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    classifiers=[c for c in CLASSIFIERS.split('\n') if c],
    keywords='Mesh Generation',
    cmake_args=cmake_args,
    setup_requires=["pytest-runner"],
    install_requires=["numpy>=1.8.0"],
    extras_require={
        "docs" : docs_require,
        "test" : tests_require,
        "build": build_require,
        "dev"  : docs_require + tests_require + build_require
    },
    zip_safe=False,
)