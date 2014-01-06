#!/usr/bin/python
# Filename: visit3D.py

version = 'v1.3'

from visit import *
import glob

MAPSTUFF_FILE = "/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/oase-mapdata.nc"
TOPO_COLORMAP = "topography"
TOPO_COLORMAP_INVERT = 0

# This module bundles python routines for handling Visit3D
# plotting more comfortably from the various visualization
# routines

# Adds clusters with names "*_infix_*.vtk" to the
# current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
#
def add_clusters(basename,infix,col_tables):
    
    # now the clusters
    cluster_pattern=basename+"*"+infix+"*.vt*"
    print "Looking for cluster files at " + cluster_pattern
    cluster_list=glob.glob(cluster_pattern)
    
    print "Processing clusters:"
    print cluster_list
    
    for cluster_file in cluster_list:
        
        # figure out cluster number for choice of color table
        # extract the number from the cluster filename
        # to select the color table based on the ID
        
        cluster_num=0
        
        try:
            #print "Extracting cluster number from " + cluster_file;
            start = cluster_file.index("_cluster_",0) + len("_cluster_")
            end = cluster_file.index(".v",start)
            cluster_num = int( cluster_file[start:end] );
            print "Plotting cluster #" + str(cluster_num)
        except ValueError:
            print "Illegal filename " + fname
            continue
        
        OpenDatabase(cluster_file)
        AddPlot("Pseudocolor", "cell_color")
    
        cp=PseudocolorAttributes();
        cp.pointSizePixels=2
        cp.legendFlag=0
        cp.lightingFlag=1
        cp.invertColorTable=0
        cp.minFlag,cp.maxFlag = 1,1
        cp.min,cp.max = 0.1, 5
        cp.opacity = 0.33
        index = cluster_num % len(col_tables);
        cp.colorTableName = col_tables[index];
        SetPlotOptions(cp)
        
    return

# Adds clusters with names "*_infix_*.vtk" to the
# current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
#
def add_clusters_with_colortable(basename,infix,color_table_name,color_count):
    
    # now the clusters
    cluster_pattern="*"+infix+"*.vt*"
    print "Looking for cluster files at " + cluster_pattern
    cluster_list=glob.glob(cluster_pattern)
    
    for cluster_file in cluster_list:
        
        # figure out cluster number for choice of color table
        # extract the number from the cluster filename
        # to select the color table based on the ID
        
        cluster_num=0
        
        try:
            start = cluster_file.index(infix,0) + len(infix)
            end = cluster_file.index(".v",start)
            cluster_num = int( cluster_file[start:end] );
        except ValueError:
            print "Illegal filename " + cluster_file
            continue
        
        OpenDatabase(cluster_file)
        AddPlot("Pseudocolor", "point_color")
        
        cp=PseudocolorAttributes();
        cp.pointSizePixels=2
        cp.legendFlag=0
        cp.lightingFlag=0
        cp.invertColorTable=0
        cp.minFlag,cp.maxFlag = 1,1
        cp.min,cp.max = 0,color_count
        cp.colorTableName = color_table_name
        cp.opacity=1
        SetPlotOptions(cp)

    return

# Adds boundaries with to the current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
# @param number of colors in the color table
#
def add_boundaries(basename,color_table_name,color_count):
    
    infix="_boundary_"
    cluster_pattern=basename+"*"+infix+"*.vt*"
    print "Looking for boundary files at " + cluster_pattern
    cluster_list=glob.glob(cluster_pattern)
    for cluster_file in cluster_list:
        cluster_num=0
        try:
            start = cluster_file.index(infix,0) + len(infix)
            end = cluster_file.index(".v",start)
            cluster_num = int( cluster_file[start:end] );
        except ValueError:
            print "Illegal filename " + cluster_file
            continue
        
        OpenDatabase(cluster_file)
        AddPlot("Pseudocolor", "point_color")
        cp=PseudocolorAttributes();
        cp.pointSizePixels = 2
        cp.legendFlag = 0
        cp.lightingFlag = 0
        cp.invertColorTable = 0
        cp.minFlag,cp.maxFlag = 1,1
        cp.min,cp.max = 0,(color_count-1)
        cp.colorTableName = color_table_name
        cp.opacity = 0.5
        SetPlotOptions(cp)
    
    return

# Add the data of one variable from the given file to
# the current window.
# @param filename
# @param variable name
# @param color table
#
def add_pseudocolor(vtk_file,variable,color_table_name,opacity,legendFlag):
    
    # open the file and add the plot
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", variable)
    p = PseudocolorAttributes()
    p.pointSizePixels = 2
    p.colorTableName = color_table_name
    p.lightingFlag = 0
    p.legendFlag = legendFlag;
    p.opacity = opacity
    SetPlotOptions(p)
    return

# Add 3D local topography
# @param "local" or "national"
#
def add_mapstuff(extent):
    
    # open the file and add the plot
    OpenDatabase(MAPSTUFF_FILE)

    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return

    # Topography
    AddPlot("Pseudocolor", extent+"_topo_3D")
    p = PseudocolorAttributes()
    p.pointSizePixels = 2
    p.colorTableName = TOPO_COLORMAP
    p.invertColorTable = TOPO_COLORMAP_INVERT
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.opacity = 1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = -100.0, 3200.0
    SetPlotOptions(p)
    
    # Rivers & Boundaries
    
    AddPlot("Pseudocolor", "as_zonal/"+extent+"_boundaries_3D")
    p = PseudocolorAttributes()
    p.colorTableName = "Greys"
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.opacity = 1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = 0, 1
    SetPlotOptions(p)
    
    AddPlot("Pseudocolor", "as_zonal/"+extent+"_rivers_3D")
    p = PseudocolorAttributes()
    p.colorTableName = "Blues"
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.invertColorTable = 1;
    p.opacity = 1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = 0, 1
    SetPlotOptions(p)
    
    return

# Closes databases connected with topo
def close_mapstuff():
    CloseDatabase(MAPSTUFF_FILE);
    return

#
# Sets default 2D view params for RADOLAN grid
# @param perspective 1,2
#
def set_view_to_radolan(perspective,scale_factor_z):

    v = GetView3D();
    
    if perspective==1:
        v.viewNormal = (0.204365, -0.63669, 0.743546)
        v.focus = (-239.212, -4222.9, 7.31354)
        v.viewUp = (-0.201314, 0.716005, 0.668438)
        v.viewAngle = 30
        v.parallelScale = 173.531
        v.nearPlane = -347.062
        v.farPlane = 347.062
        v.imagePan = (-0.00977129, 0.0399963)
        v.imageZoom = 1.4641
        v.perspective = 1
        v.eyeAngle = 2
        v.centerOfRotationSet = 0
        v.centerOfRotation = (0, 0, 0)
        v.axis3DScaleFlag = 0
        v.axis3DScales = (1, 1, 1)
        v.shear = (0, 0, 1)
    
    elif perspective==2:
        v.viewNormal = (0.996852, 0.0147052, 0.0779102)
        v.focus = (-239.212, -4222.9, 21.9406)
        v.viewUp = (-0.0779426, 0.00164413, 0.996956)
        v.viewAngle = 30
        v.parallelScale = 175.151
        v.nearPlane = -350.302
        v.farPlane = 350.302
        v.imagePan = (-0.00977129, 0.0399963)
        v.imageZoom = 1.7715
        v.perspective = 1
        v.eyeAngle = 2
        v.centerOfRotationSet = 0
        v.centerOfRotation = (0, 0, 0)
        v.axis3DScaleFlag = 0
        v.axis3DScales = (1, 1, 1)
        v.shear = (0, 0, 1)
    
    if scale_factor_z != 1.0:
        v.axis3DScaleFlag = 1
        v.axis3DScales = (1, 1, scale_factor_z)


    print "3D View Settings:"
    print v
    SetView3D(v);
    return

# Sets up standard values for axis etc
#
def set_annotations():

    a = GetAnnotationAttributes()
    a.axes3D.visible=1
    a.axes3D.autoSetScaling=0
    a.userInfoFlag=0
    a.timeInfoFlag=0
    a.legendInfoFlag=1
    a.databaseInfoFlag=1

    a.axes3D.xAxis.title.visible=0
    a.axes3D.xAxis.title.userTitle = 1
    a.axes3D.xAxis.title.userUnits = 1
    a.axes3D.xAxis.title.title = "x"
    a.axes3D.xAxis.title.units = "km"

    a.axes3D.yAxis.title.visible=0
    a.axes3D.yAxis.title.userTitle = 1
    a.axes3D.yAxis.title.userUnits = 1
    a.axes3D.yAxis.title.title = "y"
    a.axes3D.yAxis.title.units = "km"

    a.axes3D.zAxis.title.visible=0
    a.axes3D.zAxis.title.userTitle = 1
    a.axes3D.zAxis.title.userUnits = 1
    a.axes3D.zAxis.title.title = "h"
    a.axes3D.zAxis.title.units = "km"

    SetAnnotationAttributes(a)

# End of visit3D.py