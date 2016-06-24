#!BPY
"""
Name: 'GAMer mesh improvments (upy)'
Blender: 249b
Group: 'Mesh'
Tooltip: 'GAMer mesh improvments, mesh annotations and volumetric mesh generation'
"""
######################################################
# GAMer plugin using pyubic
#
# Getting GAMer:
# 
#
#
# Getting pyubic:
#
#
#
# This plugin is protected by the GPL: Gnu Public Licence
# GPL - http://www.gnu.org/copyleft/gpl.html
# Script Copyright (C) Johan Hake <hake.dev@gmail.com>
#                      Ludovic Autin <autin@scripps.edu>

# Import upy and set the UI Baseclass
import upy
upy.setUIClass()

# Importing gamer gui
from gamer.upy_gui import GamerUI

if upy.uiadaptor.host == "tk":
    from DejaVu import Viewer
    vi = Viewer()
    
    # Require a master   
    gamergui = GamerUI(master=vi, title="GAMer")
else :
    gamergui = GamerUI(title="GAMer")

# Call it
gamergui.setup()
gamergui.display()
