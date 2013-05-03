#!/usr/bin/python
import glob
import sys

# plot original data
source_file = "../../data/06.08.2011/netcdf/raa01-rx_10000-1108061550-dwd---bin.nc"
OpenDatabase( source_file )

#a = GetAnnotationAttributes()
#a.axes2D.visible=0
#SetAnnotationAttributes(a)

# Add the source data plot
AddPlot("Pseudocolor", "zh")
p = PseudocolorAttributes()
p.colorTableName="hot_desaturated"
p.invertColorTable=0
p.legendFlag=1
p.minFlag, p.maxFlag = 1,1
p.min,p.max = 0,50
p.opacity=0.33
SetPlotOptions(p)

# Add threshold operator
AddOperator("Threshold")
t = ThresholdAttributes();
t.lowerBounds=(10)
print t
#t.lbound = 10
SetOperatorOptions(t)

DrawPlots();

# plot the clusters
cluster_basename = "Release/*_cluster_*.vtk"
cluster_basename = "Debug/*_cluster_*.vtk"
list = glob.glob( cluster_basename )

count = 0;
#col_tables = ["Purples","Blues","Oranges","Greens","Reds","Set1","Set2","Set3","hot","hot_and_cold","hot_desaturated"]
col_tables = ["Purples","Blues","Oranges","Reds","RdYlGn"]
for fname in list:

    # add plot                                                                                                                                                                      
    OpenDatabase(fname);
    AddPlot("Pseudocolor", "cluster")

    # set plot attributes accordingly                                                                                                                                               
    cp=PseudocolorAttributes();
    cp.minFlag,cp.maxFlag = 1,1
    cp.min,cp.max = 0.0, 50.0
    cp.pointSizePixels=10
    cp.legendFlag, cp.lightingFlag = 0,0
    cp.invertColorTable=0
    index = count % len(col_tables);
    cp.colorTableName = col_tables[index];
    cp.opacity=0.125

    SetPlotOptions(cp)
    count = count + 1;

DrawPlots()
