# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2021
# by Christopher T. Lee and contributors
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http:#www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# ***************************************************************************

import sys

bver_to_pyver = {
    # "2.79": '3.5',
    "2.83": '3.7',
    "2.93": '3.9',
    "3.3": "3.10",
    "3.6": "3.10",

}

if __name__ == "__main__":
    print(bver_to_pyver[sys.argv[1]])
