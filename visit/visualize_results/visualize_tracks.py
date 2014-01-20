#!/usr/bin/python

# Placeholders for substitution

MEANIE3D_HOME     = "P_M3D_HOME"
SOURCE_DIR        = "P_SOURCE_DIR"
BASENAME          = "P_BASENAME"

# Appending the module path is crucial

sys.path.append(MEANIE3D_HOME+"/visit/modules")

import glob
import visit2D
import visitUtils

# Modify view parameters
visit2D.set_view_to_radolan();

# TODO: plot the underlying map data

# Plot the Tracks

if BASENAME == "":
    track_pattern = SOURCE_DIR + "/*-track_*.vtk"
else:
    track_pattern = SOURCE_DIR + "/" + BASENAME + "-track_*.vtk"

list = glob.glob(track_pattern)

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
