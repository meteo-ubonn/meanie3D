#!/usr/bin/python
import glob
import sys

# plot original data
source_file = "../Debug/raa01-rx_10000-1106221500-dwd---bin.vtk"
OpenDatabase( source_file )

a = AnnotationAttributes()
print a

AddPlot("Pseudocolor", "reflectivity")
p = PseudocolorAttributes()
p.colorTableName = "xray"
p.legendFlag=0
p.lightingFlag=0
p.invertColorTable=0
#p.minFlag,p.maxFlag = 1,1
#p.min,p.max = 0.0, 50.0
SetPlotOptions(p)
DrawPlots();

# plot the clusters
#cluster_basename = "Release/*_cluster_*.vtk"
cluster_basename = "../Debug/*_cluster_*.vtk"
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
    cp.pointSizePixels=8
    cp.legendFlag, cp.lightingFlag = 0,0
    cp.invertColorTable=0
    index = count % len(col_tables);
    cp.colorTableName = col_tables[index];
    cp.opacity=0.02

    SetPlotOptions(cp)
    count = count + 1;

DrawPlots()