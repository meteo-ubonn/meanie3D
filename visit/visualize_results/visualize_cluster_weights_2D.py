#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME = "M3D_HOME_P"
SOURCE_FILE   = "SOURCE_FILE_P"
VAR_NAME      = "VAR_NAME_P"

# Appending the module path is crucial

sys.path.append(MEANIE3D_HOME+"/visit/modules")

import glob
import os
import sys
import time
import visit2D
import visitUtils
from subprocess import call

# plot original data

# Modify view parameters                                                                                                                                                                                         
visit2D.set_view_to_radolan()

# plot data set with opacity

OpenDatabase(SOURCE_FILE)

AddPlot("Pseudocolor", VAR_NAME)
p = PseudocolorAttributes()
p.colorTableName = "hot_desaturated"
p.legendFlag=1
p.lightingFlag=1
p.invertColorTable=0
p.pointSizePixels=4
SetPlotOptions(p)
DrawPlots();

# plot weights?

if "no" == "yes":
    weights_file="../Debug/*-weights.vtk"
    list = glob.glob(weights_file)
    OpenDatabase(list[0])

    AddPlot("Pseudocolor", "weight")
    p = PseudocolorAttributes()
    p.colorTableName = "hot_desaturated"
    p.legendFlag=1
    p.lightingFlag=0
    p.invertColorTable=0
    p.pointSizePixels=3
    SetPlotOptions(p)
    DrawPlots();

# plot the clusters
#cluster_basename = "Release/*_cluster_*.vtk"
cluster_basename = "../Debug/*-clusters_weight*.vtk"
list = glob.glob( cluster_basename )

count = 0;
for fname in list:

    # add plot                                                                                                                                                                      
    OpenDatabase(fname);
    AddPlot("Pseudocolor", "weight")

    # set plot attributes accordingly                                                                                                                                               
    cp=PseudocolorAttributes();
    cp.pointSizePixels=3
    cp.opacity=0.5
    cp.minFlag=1
    cp.maxFlag=1
    cp.min=0
    cp.max=0.5
    cp.legendFlag=0
    cp.colorTableName = "rainbow";
    if count%2==0:
        cp.invertColorTable=1
    else:
        cp.invertColorTable=0
    SetPlotOptions(cp)

    count = count+1

# plot vectors

OpenDatabase("../Debug/meanshift_vectors-spatial.vtk")
AddPlot("Vector","vectors")
vp = VectorAttributes()
vp.useStride = 1
vp.stride = 1
vp.scale = 1
vp.scaleByMagnitude = 1
vp.autoScale = 0
vp.headSize = 0.33
vp.headOn = 1
vp.colorByMag = 1
vp.useLegend = 1
vp.stemWidth = 0.08
vp.origOnly = 1
vp.colorTableName = "gray"
vp.invertColorTable = 1
SetPlotOptions(vp)

DrawPlots()

quit()