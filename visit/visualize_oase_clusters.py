#!/usr/bin/python
import glob
import sys
import os

# plot original data
source_file = "SOURCE_FILE"

fileName, fileExtension = os.path.splitext(source_file)

print "Plotting " + source_file
OpenDatabase( source_file )

#a = GetAnnotationAttributes()
#a.axes2D.visible=0
#SetAnnotationAttributes(a)

# Modify view parameters
v = GetView2D()
v.windowCoords = (-418.462, 292.538, -4446.64, -3759.64)
v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
SetView2D(v)

# Add the source data plot
#AddPlot("Pseudocolor", "cband_radolan_rx")
AddPlot("Pseudocolor", "zh")
p = PseudocolorAttributes()
p.colorTableName="xray"
p.invertColorTable=1
p.pointSizePixels=5
p.legendFlag=1
p.minFlag, p.maxFlag = 1,1
p.min,p.max = 0,60
p.opacity=1.0
SetPlotOptions(p)

# Add threshold operator
AddOperator("Threshold")
t = ThresholdAttributes();
t.lowerBounds=(10)
print t
#t.lbound = 10
SetOperatorOptions(t)

#AddOperator("ThreeSlice");
#s=ThreeSliceAttributes();
#s.x = -231.9620
#s.y = -4180.1450
#s.z = 0.0
#SetOperatorOptions(s);

DrawPlots();

# plot the clusters
# cluster_basename = "Release/*_cluster_*.vtk"
list = glob.glob( "../Debug/*_cluster_*.vtk" )

count = 0;
col_tables = ["Purples","Blues","Oranges","Greens","Reds","Paired"]
#col_tables = ["Purples","Blues","Oranges","Reds","RdYlGn","Orangehot"]
for fname in list:

    # add plot                                                                                                                                                                      
    OpenDatabase(fname);
    AddPlot("Pseudocolor", "cluster")

    # set plot attributes accordingly                                                                                                                                               
    cp=PseudocolorAttributes();
    #cp.minFlag,cp.maxFlag = 1,1
    #cp.min,cp.max = 0.0, 50.0
    cp.pointSizePixels=2
    cp.legendFlag, cp.lightingFlag = 0,0
    cp.invertColorTable=0
    index = count % len(col_tables);
    cp.colorTableName = col_tables[index];
    cp.opacity=0.5

    SetPlotOptions(cp)
    count = count + 1;

list = glob.glob( "../Debug/*clusters_centers.vtk" )
OpenDatabase(list[0])
AddPlot("Label","geometrical_center")
a = LabelAttributes()
#a.legendFlag = 1
#a.showNodes = 0
#a.showCells = 1
a.restrictNumberOfLabels = 0
#a.numberOfLabels = 200
a.specifyTextColor1 = 0
a.textColor1 = (255, 0, 0, 0)
a.textHeight1 = 0.03
a.legendFlag = 0
#a.specifyTextColor2 = 0
#a.textColor2 = (0, 0, 255, 0)
#a.textHeight2 = 0.03
a.formatTemplate = "%g"
SetPlotOptions(a)


DrawPlots()
