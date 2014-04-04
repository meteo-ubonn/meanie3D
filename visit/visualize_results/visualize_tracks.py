#!/usr/bin/python

# Placeholders for substitution

MEANIE3D_HOME     = "P_M3D_HOME"
SOURCE_DIR        = "P_SOURCE_DIR"
BASENAME          = "P_BASENAME"

# Appending the module path is crucial

sys.path.append(MEANIE3D_HOME+"/visit/modules")

# Import modules
import sys
sys.path.append(MEANIE3D_HOME+"/visit/modules")
import glob
import os
import time
import visit2D
import visitUtils
from subprocess import call

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes

a = GetAnnotationAttributes()
a.axes2D.visible=1
a.axes2D.autoSetScaling=0
a.axes2D.xAxis.title.visible=0
a.axes2D.yAxis.title.visible=0
a.legendInfoFlag=1
a.databaseInfoFlag=0
a.userInfoFlag=0
a.timeInfoFlag=0
SetAnnotationAttributes(a)

# Add gray/black background gradient
visitUtils.add_background_gradient();

print "Cleaning up *.vtk *.vtr *.png"
return_code=call("rm -f *.vtk *.vtr *.png", shell=True)

# Set view to nationwide composite
visit2D.set_view_to_radolan();

print "-- Creating colortables ---"
num_colors = visitUtils.create_cluster_colortable("cluster_colors")
visitUtils.create_topography_colortable()
print "    done."

# Plot the Tracks

if BASENAME == "":
    track_pattern = SOURCE_DIR + "/*-track_*.vtk"
else:
    track_pattern = SOURCE_DIR + "/" + BASENAME + "-track_*.vtk"

# add topograpy
visit2D.add_topography("national_topo_2D")

list = sorted(glob.glob(track_pattern))
count = 0;
for fname in list:

    # add plot                                                                                                                                                                      
    OpenDatabase(fname);
    AddPlot("Pseudocolor", "track_step")

    # set plot attributes accordingly                                                                                                                                               
    cp=PseudocolorAttributes();
    cp.pointSizePixels=10
    
    # min/max could be subject to parametrization
    cp.minFlag,cp.maxFlag=1,1
    cp.min,cp.max=0,50
    
    # plot the legend for the first one
    if count==0:
        cp.legendFlag=1
    else:
        cp.legendFlag=0

    cp.colorTableName = "hot";

    SetPlotOptions(cp)

    count = count+1

DrawPlots()
