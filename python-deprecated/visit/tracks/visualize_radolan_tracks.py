#!/usr/bin/python

# Placeholders for substitution

MEANIE3D_HOME     = "P_M3D_HOME"
SOURCE_DIR        = "P_SOURCE_DIR"
BASENAME          = "P_BASENAME"

# Import modules
import sys
sys.path.append(MEANIE3D_HOME+"/python/python-modules")
import glob
from meanie3D import visit2D
from meanie3D.app import utils
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
# visitUtils.add_background_gradient();

# Set view to nationwide composite
visit2D.set_view_to_radolan();

#print "-- Creating colortables ---"
#num_colors = visitUtils.create_cluster_colortable("cluster_colors")
#visitUtils.create_topography_colortable()
#print "    done."

# Plot the Tracks

if BASENAME == "":
    track_pattern = SOURCE_DIR + "/*-track_*.vtk"
else:
    track_pattern = SOURCE_DIR + "/" + BASENAME + "-track_*.vtk"

# add topograpy
#visit2D.add_topography("national_topo_2D")

# add 2D topograpy
visit2D.add_map_borders("national")
visit2D.add_map_rivers("national")

list = sorted(glob.glob(track_pattern))

print "Looking with pattern " + track_pattern
print "Found the following track files:"
print list

count = 0;
for fname in list:

    # add plot                                                                                                                                                                      
    OpenDatabase(fname);
    AddPlot("Pseudocolor", "track_step")

    # set plot attributes accordingly                                                                                                                                               
    cp=PseudocolorAttributes();
    cp.pointSizePixels=5
    
    # min/max could be subject to parametrization
    cp.minFlag,cp.maxFlag=1,1
    cp.min,cp.max=0,50
    
    cp.legendFlag=0
    
    # plot the legend for the first one
    #if count==0:
    #    cp.legendFlag=1
    #else:
    #    cp.legendFlag=0

    cp.colorTableName = "hot_desaturated";

    SetPlotOptions(cp)

    count = count+1

DrawPlots()


utils.save_window("tracks",0)

print "Cleaning up *.vtk"
return_code=call("rm -f *.vtk", shell=True)

quit()
