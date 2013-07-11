#!/usr/bin/python
import glob
import sys

# plot original data
source_file = "SOURCE_FILE"

weights_file="../Debug/*-weights.vtk"
list = glob.glob(weights_file)
OpenDatabase(list[0])

AddPlot("Pseudocolor", "weight")
p = PseudocolorAttributes()
p.colorTableName = "hot_desaturated"
p.legendFlag=1
p.lightingFlag=0
p.invertColorTable=0
p.pointSizePixels=20
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
    cp.minFlag,cp.maxFlag = 1,1
    #cp.min,cp.max = 0.0, 50.0
    cp.pointSizePixels=5
    if count == 0:
        cp.legendFlag, cp.lightingFlag = 1,0
    else:
        cp.legendFlag, cp.lightingFlag = 0,0

    cp.invertColorTable=0
    cp.colorTableName = "hot_desaturated";
    #cp.opacity=0.02

    SetPlotOptions(cp)
    count = count + 1;

    #if count > 200:
    #    break

OpenDatabase("../Debug/meanshift_vectors-spatial.vtk")
AddPlot("Vector","vectors")
vp = VectorAttributes()
vp.useStride = 1
vp.stride = 1
vp.scale = 1
vp.scaleByMagnitude = 1
vp.autoScale = 0
vp.headSize = 0.25
vp.headOn = 1
vp.colorByMag = 0
vp.useLegend = 1
vp.stemWidth = 0.08
vp.origOnly = 1
SetPlotOptions(vp)

DrawPlots()
