#!/usr/bin/python
# Filename: visit3D.py

version = '0.1'

from visit import *
import glob

TOPO_FILE = "/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/oase-mapdata.nc"
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

# Add 3D topography
#
def add_topography(var_name):
    
    # open the file and add the plot
    OpenDatabase(TOPO_FILE)
    AddPlot("Pseudocolor", var_name)
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
    return

# Closes databases connected with topo
def close_topography():
    CloseDatabase(TOPO_FILE);
    return

#
# Sets default 2D view params for RADOLAN grid
#
def set_view_to_radolan():

    v = GetView3D();

#    v.viewNormal = (0.656802,-0.498223,0.566025)
#    v.focus = (-239.212,-4222.9,7.375)
#    v.viewUp = (-0.457525,0.333371,0.824339)
#    v.viewAngle = 30
#    v.parallelScale = 173.528
#    v.nearPlane = -347.056
#    v.farPlane = 347.056
#    v.imagePan = (0, 0)
#    v.imageZoom = 1.4641
#    v.perspective = 1
#    v.eyeAngle = 2
#    v.centerOfRotationSet = 0
#    v.centerOfRotation = (0, 0, 0)
#    v.axis3DScaleFlag = 0
#    v.axis3DScales = (1, 1, 1)
#    v.shear = (0, 0, 1)

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

    print "3D View Settings:"
    print v
    SetView3D(v);
    return


# End of visit3D.py