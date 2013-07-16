#!/usr/bin/python
import glob
import sys

# plot original data
source_file = "../Debug/raa01-rx_10000-1307011600-dwd---bin_reflectivity.vtk"

# Modify view parameters                                                                                                                                                                                         
v = GetView2D()
v.windowCoords = (-244.462, 376.538, -4551.64, -3806.64)
v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
SetView2D(v)

# plot data set with opacity

OpenDatabase(source_file)

AddPlot("Pseudocolor", "reflectivity")
p = PseudocolorAttributes()
p.colorTableName = "hot_desaturated"
p.legendFlag=1
p.lightingFlag=0
p.invertColorTable=0
p.pointSizePixels=5
p.opacity=0.25
SetPlotOptions(p)
DrawPlots();

# plot weights

weights_file="../Debug/*-weights.vtk"
list = glob.glob(weights_file)
#OpenDatabase(list[0])

#AddPlot("Pseudocolor", "weight")
#p = PseudocolorAttributes()
#p.colorTableName = "hot_desaturated"
#p.legendFlag=1
#p.lightingFlag=0
#p.invertColorTable=0
#p.pointSizePixels=20
#SetPlotOptions(p)
#DrawPlots();


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
    cp.pointSizePixels=8
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
vp.headSize = 0.1
vp.headOn = 1
vp.colorByMag = 1
vp.useLegend = 1
vp.stemWidth = 0.08
vp.origOnly = 1
vp.colorTableName = "gray"
vp.invertColorTable = 1
SetPlotOptions(vp)

DrawPlots()
