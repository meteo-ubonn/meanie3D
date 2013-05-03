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
def add_clusters_2D(basename,infix,col_tables):

    cluster_pattern=basename+"*" + infix + "*.vtk"
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
            print "Extracting cluster number from " + cluster_file;
            start = cluster_file.index(infix,0) + len(infix)
            end = cluster_file.index(".vtk",start)
            cluster_num = int( cluster_file[start:end] );
            print "Plotting cluster #" + str(cluster_num)
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


# Closes the databases associated with the cluster files
# parameters follow the same rules as add_clusters_2D
def close_clusters_2D(basename,infix):
    cluster_pattern=basename+"*" + infix + "_" + "*.vtk"
    cluster_list=glob.glob(cluster_pattern)
    for cluster_file in cluster_list:
        CloseDatabase(cluster_file)
    return


# Add the data of one variable from the given file to
# the current window.
# @param filename
# @param variable name
# @param color table
#
def add_pseudocolor_2D(vtk_file,variable,color_table_name):

    # open the file and add the plot
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", variable)
    
    p = PseudocolorAttributes()
    p.colorTableName = color_table_name
    #p.legendFlag=1
    p.lightingFlag=1
    #p.invertColorTable=1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = 0.0, 50.0
    p.opacity=1.0
    SetPlotOptions(p)
    return

def add_modes(modes_file):
    OpenDatabase(modes_file)
    AddPlot("Label","mode")
    a = LabelAttributes()
    #a.legendFlag = 1
    #a.showNodes = 0
    #a.showCells = 1
    a.restrictNumberOfLabels = 0
    #a.numberOfLabels = 200
    a.specifyTextColor1 = 0
    a.textColor1 = (255, 0, 0, 0)
    a.textHeight1 = 0.03
    #a.specifyTextColor2 = 0
    #a.textColor2 = (0, 0, 255, 0)
    #a.textHeight2 = 0.03
    a.formatTemplate = "%g"
    SetPlotOptions(a)
    return

# Saves PNG file
# @param basename for the file
# @param is numbering appended?
def save_window(basename,progressive):
    s = GetSaveWindowAttributes()
    s.outputToCurrentDirectory = 1
    s.outputDirectory = "."
    s.fileName=basename
    s.width = 1024
    s.height = 1024
    s.quality = 100
    s.progressive=progressive
    SetSaveWindowAttributes(s)
    SaveWindow()
    return

# End of visit2D.py