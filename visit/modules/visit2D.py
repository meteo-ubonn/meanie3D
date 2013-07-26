#!/usr/bin/python
# Filename: visit2D.py

version = '0.1'

from visit import *
import glob

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

    cluster_pattern=basename+"*" + infix + "*.vtk"
    #print "Looking for cluster files at " + cluster_pattern
    cluster_list=glob.glob(cluster_pattern)
    
    #print "Processing clusters:"
    #print cluster_list
    
    for cluster_file in cluster_list:
        
        # figure out cluster number for choice of color table
        # extract the number from the cluster filename
        # to select the color table based on the ID
        
        cluster_num=0
        
        try:
            #print "Extracting cluster number from " + cluster_file;
            start = cluster_file.index(infix,0) + len(infix)
            end = cluster_file.index(".vtk",start)
            cluster_num = int( cluster_file[start:end] );
            #print "Plotting cluster #" + str(cluster_num)
        except ValueError:
            print "Illegal filename " + fname
            continue
        
        OpenDatabase(cluster_file)
        AddPlot("Pseudocolor", "cluster")

        cp=PseudocolorAttributes();
        #cp.minFlag,cp.maxFlag = 1,1
        #cp.min,cp.max = 10.0, 60.0
        cp.pointSizePixels=5
        cp.legendFlag=0
        cp.lightingFlag=1
        cp.invertColorTable=0
        index = cluster_num % len(col_tables);
        cp.colorTableName = col_tables[index];
        cp.opacity=0.05
        SetPlotOptions(cp)
    return


# Add the data of one variable from the given file to
# the current window.
# @param filename
# @param variable name
# @param color table
#
def add_pseudocolor(vtk_file,variable,color_table_name,lf):

    # open the file and add the plot
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", variable)
    
    p = PseudocolorAttributes()
    p.colorTableName = color_table_name
    p.legendFlag=lf
    p.lightingFlag = 1
    p.opacity=1.0
    SetPlotOptions(p)
    return

# End of visit2D.py