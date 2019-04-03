# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2018
# by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
#    and Michael Holst

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# ***************************************************************************

include(ExternalProject)
ExternalProject_Add(
    casc
    URL https://github.com/ctlee/casc/archive/development.tar.gz
    #URL_HASH SHA1=34cfe1ae39c4abd63f88e0056202d376c01ff190
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/casc
    INSTALL_COMMAND ""
)

add_subdirectory(tetgen)
# add_subdirectory(triangle)
