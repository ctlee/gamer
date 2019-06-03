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

import skbuild
from skbuild import setup
from distutils.version import LooseVersion
from setuptools import find_packages

docs_require = ["sphinx", "breathe", "exhale", "sphinx-issues"]
tests_require = ["pytest"]

setup(
    name='pygamer',
    version='2.0.1',
    author='Christopher T. Lee',
    author_email='ctlee@ucsd.edu',
    description='GAMer: Geometry-preserving Adaptive Mesher',
    long_description='Python wrapper around the GAMer C++ library for mesh generation.',
    cmake_args=['-DBUILD_PYGAMER=ON'],
    zip_safe=False,
    packages=find_packages(),
    setup_requires=["pytest-runner"],
    extras_require={
        "docs" : docs_require,
        "test" : tests_require,
        "dev"  : docs_require + tests_require
    },
)