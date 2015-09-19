#!/usr/bin/python

# ------------------------------------------------------------------------------
# Filename: visit2D.py
#
# This module bundles python routines for handling Visit2D plotting
#
# @author Juergen Simon (juergen_simon@mac.com)
# ------------------------------------------------------------------------------

from visit import *
import sys
sys.path.append(".")

import glob
import os
import time
import visitUtils
from pprint import pprint
from subprocess import call

# ------------------------------------------------------------------------------
# The following variables are important for referencing the mapstuff and
# topography data
# ------------------------------------------------------------------------------
MAPS_HOME=os.environ['MEANIE3D_MAPDATA']
TOPO_FILE = MAPS_HOME+"/oase-mapdata.nc"
TOPO_COLORMAP = "topography"
TOPO_COLORMAP_INVERT = 0

# ------------------------------------------------------------------------------
# Adds clusters with names "*_infix_*.vtk" to the
# current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
# ------------------------------------------------------------------------------
def add_clusters(basename,infix,col_tables,opacity):

    cluster_pattern=basename+"*" + infix + "*.vtk"
    #print "Looking for cluster files at " + cluster_pattern
    cluster_list=glob.glob(cluster_pattern)
    
    #print "Processing clusters:"
    #print cluster_lis
    
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
        cp.pointSizePixels=2
        cp.legendFlag=0
        cp.lightingFlag=1
        cp.invertColorTable=0
        index = cluster_num % len(col_tables);
        cp.colorTableName = col_tables[index];
        cp.opacity=opacity
        SetPlotOptions(cp)
    return

# ------------------------------------------------------------------------------
# Adds clusters with names "*_infix_*.vtk" to the
# current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
# ------------------------------------------------------------------------------
def add_clusters_with_colortable(basename,infix,color_table_name,color_count):
    
    # now the clusters
    cluster_pattern="*"+infix+"*.vt*"
    print "Looking for cluster files at " + cluster_pattern
    cluster_list=sorted(glob.glob(cluster_pattern))
    
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
        
        SetActivePlots(GetNumPlots()-1)
        SetPlotOptions(cp)
    
    return

# ------------------------------------------------------------------------------
# Add the data of one variable from the given file to
# the current window.
# @param filename
# @param variable name
# @param color table
# ------------------------------------------------------------------------------
def add_pseudocolor(vtk_file,variable,color_table_name,opacity,legendFlag):

    # open the file and add the plot
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", variable)
    p = PseudocolorAttributes()
    p.colorTableName = color_table_name
    p.legendFlag=legendFlag
    p.lightingFlag=1
    p.opacity=opacity
    # p.minFlag,p.maxFlag = 1,1
    # p.min,p.max = 35,92.5
    SetPlotOptions(p)
    return

# ------------------------------------------------------------------------------
# Sets default 2D view params for RADOLAN grid
# ------------------------------------------------------------------------------
def set_view_to_radolan():
    v = GetView2D()
    v.windowCoords = (-523.962, 376.038, -4659.15, -3759.15)
    v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    SetView2D(v)
    return

# ------------------------------------------------------------------------------
# Sets 2D view params 
# @param dictionary with 'windowCoords' and 'viewportCoords' keys
# ------------------------------------------------------------------------------
def set_view(params):
    v = GetView2D()
    v.viewportCoords = tuple(params["viewportCoords"])
    v.windowCoords = tuple(params["windowCoords"])
    SetView2D(v)
    return

# ------------------------------------------------------------------------------
# Add 3D topography
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Closes databases connected with topo
# ------------------------------------------------------------------------------
def close_topography():
    CloseDatabase(TOPO_FILE);
    return

# ------------------------------------------------------------------------------
# Add 2D local topography, rivers and borders
# @param "local" or "national"
# ------------------------------------------------------------------------------
def add_mapstuff(extent):
    
    # Topography
    add_map_topography(extent);
    
    # Rivers & Boundaries
    add_map_borders(extent);
    add_map_rivers(extent);
    
    return

# ------------------------------------------------------------------------------
# Add 2D local topography
# @param "local" or "national"
# ------------------------------------------------------------------------------
def add_map_topography(extent):
    
    # open the file and add the plot
    OpenDatabase(TOPO_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
    # Topography
    AddPlot("Pseudocolor", extent+"_topo_2D")
    p = PseudocolorAttributes()
    p.pointSizePixels = 2
    #p.colorTableName = TOPO_COLORMAP
    #p.invertColorTable = TOPO_COLORMAP_INVERT
    p.invertColorTable=1
    p.colorTableName = "xray"
    p.scaling = p.Skew
    p.skewFactor=0.0001
    
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.opacity = 1
    #p.minFlag,p.maxFlag = 1,1
    #p.min,p.max = -100.0, 3200.0
    SetPlotOptions(p)

    return

# ------------------------------------------------------------------------------
# Add 2D borders
# @param "local" or "national"
# ------------------------------------------------------------------------------
def add_map_borders(extent):
    
    # open the file and add the plot
    OpenDatabase(TOPO_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
    # Rivers & Boundaries
    
    AddPlot("Pseudocolor", "as_zonal/"+extent+"_boundaries_2D")
    p = PseudocolorAttributes()
    p.colorTableName = "Greys"
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.opacity = 1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = 0, 1
    SetPlotOptions(p)
    
    return

# ------------------------------------------------------------------------------
# Add 2D borders
# @param "local" or "national"
# ------------------------------------------------------------------------------
def add_map_rivers(extent):
    
    # open the file and add the plot
    OpenDatabase(TOPO_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
    # Rivers & Boundaries
    
    AddPlot("Pseudocolor", "as_zonal/"+extent+"_rivers_2D")
    p = PseudocolorAttributes()
    p.colorTableName = "hot"
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.invertColorTable = 1;
    p.opacity = 1
    p.minFlag,p.maxFlag = 1,1
    p.min,p.max = 0, 1
    SetPlotOptions(p)
    
    return

# ------------------------------------------------------------------------------
# Sets up standard values for axis etc
# ------------------------------------------------------------------------------
def set_annotations():
    
    a = GetAnnotationAttributes()
    a.axes2D.visible=1
    a.axesArray.autoSetScaling=0
    a.axes2D.xAxis.title.visible=0
    a.axes2D.yAxis.title.visible=0
    a.legendInfoFlag=1
    a.userInfoFlag=0
    a.timeInfoFlag=0
    a.databaseInfoFlag=1

    a.axes2D.xAxis.title.visible=0
    a.axes2D.xAxis.title.userTitle = 1
    a.axes2D.xAxis.title.userUnits = 1
    a.axes2D.xAxis.title.title = "x"
    a.axes2D.xAxis.title.units = "km"
    
    a.axes2D.yAxis.title.visible=0
    a.axes2D.yAxis.title.userTitle = 1
    a.axes2D.yAxis.title.userUnits = 1
    a.axes2D.yAxis.title.title = "y"
    a.axes2D.yAxis.title.units = "km"
    
    SetAnnotationAttributes(a)

# ------------------------------------------------------------------------------
# Generic routine for visualizing 3D clusters in two perspectives
#
# TODO: control perspectives via configuration options
#
# The following configuration options exist:
#
# 'NETCDF_DIR' : directory with the source data files
# 'CLUSTER_DIR' : directory with the cluster results
# 'M3D_HOME' : home directory of meanie3D (for the mapstuff file and modules)
# 'RESUME' : if true, the existing image files are not wiped and work is
#            picked up where it left off. Otherwise all existing images
#            are deleted and things are started from scratch
# 'WITH_BACKGROUND_GRADIENT' : add a gray background gradient to the canvas?
# 'WITH_TOPOGRAPHY' : use the topography data from the mapstuff file?
# 'WITH_RIVERS_AND_BOUNDARIES' : add rivers and boundaries?
# 'WITH_SOURCE_BACKROUND' : re-add the source data when plotting clusters?
# 'WITH_DATETIME' : add a date/time label?
# 'CREATE_SOURCE_MOVIE' : create a movie from the source images?
# 'CREATE_CLUSTERS_MOVIE' : create a movie from the cluster images?
# 'GRID_EXTENT' : "national" or "local"
# 'CONVERSION_PARAMS' : parameters for meanie3D-cfm2vtk
# 'VARIABLES' : list of variables for the source data
# 'LOWER_TRESHOLDS' : bottom cutoff for each variable,
# 'UPPER_TRESHOLDS' : top cutoff for each variable,
# 'VAR_MIN' : lowest value on legend
# 'VAR_MAX' : highest value on legend
# 'COLORTABLES' : colortable to use for each variable
# 'COLORTABLES_INVERT_FLAGS' : flag indicating inversion of colortable
#                              for each variable
# 'OPACITY' : opacity to use for each variable
# ------------------------------------------------------------------------------
def visualization(conf):

    print "-------------------------------------------------"
    print "visit2D.visualization started with configuration:"
    print "-------------------------------------------------"
    pprint(conf)
    print "-------------------------------------------------"

    bin_prefix = "export DYLD_LIBRARY_PATH="+visitUtils.get_dyld_library_path()+";"
    conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"

    # Silent
    SuppressMessages(True)
    SuppressQueryOutputOn()

    # Set view and annotation attributes

    print "Setting annotation attributes:"
    set_annotations()

    if conf['RESUME'] == False:
        print "Removing results from previous runs"
        return_code=call("rm -f images movies *.vtk *.vtr *.png", shell=True)
    else:
        print "Removing intermediary files from previous runs"
        return_code=call("rm -f *.vtk *.vtr", shell=True)

    # Set view to nationwide composite
    set_view_to_radolan();

    print "-- Creating colortables ---"
    num_colors = visitUtils.create_cluster_colortable("cluster_colors")

    if conf['WITH_TOPOGRAPHY']:
        visitUtils.create_topography_colortable()

    if conf['WITH_BACKGROUND_GRADIENT']:
        visitUtils.add_background_gradient();

    # Glob the netcdf directory
    netcdf_files = sorted(glob.glob(conf['NETCDF_DIR']+"/*.nc"));

    print "Processing files in directory " + conf['NETCDF_DIR']

    # Keep track of number of images to allow
    # forced re-set in time to circumvent the
    # Visit memory leak
    image_count=0

    # optional method of setting window for 2D
    if 'VIEW' in conf.keys():
        set_view(conf['VIEW'])
    
    for netcdf_file in netcdf_files:
        
        # construct the cluster filename and find it
        # in the cluster directory
        
        netcdf_path,filename    = os.path.split(netcdf_file);
        basename                = os.path.splitext(filename)[0]
        
        cluster_file            = conf['CLUSTER_DIR']+"/"+basename+"-clusters.nc"
        label_file              = basename+"-clusters_centers.vtk"
        
        # check if the files both exist
        print "Visualzing file "+netcdf_file+" and cluster file "+cluster_file
        if not os.path.exists(cluster_file):
            print "Cluster file does not exist. Skipping."
            continue
        
        # predict the filenames for checking on resume
        number_postfix = str(image_count).rjust(4,'0') + ".png";
        
        source_open = False
        skip_source = False
        
        if conf['RESUME'] == True:
            fn = "source_" + number_postfix
            if os.path.exists(fn):
                print "Skipping existing file " + fn
                skip_source = True
        
        if skip_source == False:

            if conf['CREATE_SOURCE_MOVIE']:
    
                OpenDatabase(netcdf_file);
                source_open = True

                if conf['WITH_TOPOGRAPHY']:
                    print "-- Adding topography data --"
                    add_topography("national_topo_2D")
                
                if conf['WITH_RIVERS_AND_BOUNDARIES']:
                    print "-- Adding map data --"
                    add_map_rivers("national")
                    add_map_borders("national")
                
                if conf['WITH_DATETIME']:
                    print "-- Adding timestamp --"
                    visitUtils.add_datetime(netcdf_file)
                
                print "-- Plotting source data --"
                start_time = time.time()
                
                variables = conf['VARIABLES']
                colortables = conf['COLORTABLES']
                opacities = conf['OPACITY']
                
                for i in range(len(variables)):

                    add_pseudocolor(netcdf_file, variables[i], colortables[i], opacities[i],1)

                    p = PseudocolorAttributes()
                    p.colorTableName = str(colortables[i])
                    p.minFlag,p.maxFlag=1,1
                    p.min,p.max=conf['VAR_MIN'][i],conf['VAR_MAX'][i]
                    p.invertColorTable = conf['COLORTABLES_INVERT_FLAGS'][i]
                    SetPlotOptions(p)
                    
                    AddOperator("Threshold")
                    t = ThresholdAttributes();
                    t.lowerBounds=(conf['LOWER_TRESHOLDS'][i])
                    t.upperBounds=(conf['UPPER_TRESHOLDS'][i])
                    SetOperatorOptions(t)

                DrawPlots();
                visitUtils.save_window("source_",1)

                DeleteAllPlots()
                ClearWindow()
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
        
        if conf['CREATE_CLUSTERS_MOVIE']:
            
            skip = False
                
            if conf['RESUME'] == True:
                fn = "tracking_" + number_postfix
                
                if os.path.exists(fn):
                    print "Skipping existing files " + fn
                    skip = True

            if skip == False:

                start_time = time.time()
                print "-- Converting clusters to .vtr --"
                
                # build the clustering command
                command=conversion_bin+" -f "+cluster_file+" "+conf['CONVERSION_PARAMS']
                print command
                return_code = call( command, shell=True)
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
                
                if conf['WITH_TOPOGRAPHY']:
                    print "-- Adding topography data --"
                    add_topography("national_topo_2D")
                
                if conf['WITH_RIVERS_AND_BOUNDARIES']:
                    print "-- Adding map data --"
                    add_map_rivers("national")
                    add_map_borders("national")
                
                if conf['WITH_DATETIME']:
                    print "-- Adding timestamp --"
                    visitUtils.add_datetime(netcdf_file)
                
                print "-- Rendering cluster scene --"
                start_time = time.time()
                
                if conf['WITH_SOURCE_BACKROUND']:
                    
                    if not source_open:
                        OpenDatabase(netcdf_file);
                        source_open = True
                    
                    for i in range(len(variables)):
                        
                        add_pseudocolor(netcdf_file,
                                        conf['VARIABLES'][i],
                                        conf['COLORTABLES'][i],
                                        conf['OPACITY'][i],1)
                                        
                        p = PseudocolorAttributes()
                        p.invertColorTable = conf['COLORTABLES_INVERT_FLAGS'][i]
                        p.minFlag,p.maxFlag=1,1
                        p.min,p.max=conf['VAR_MIN'][i],conf['VAR_MAX'][i]
                        SetPlotOptions(p)
                        
                        AddOperator("Threshold")
                        t = ThresholdAttributes();
                        t.lowerBounds=(conf['LOWER_TRESHOLDS'][i])
                        t.upperBounds=(conf['UPPER_TRESHOLDS'][i])
                        SetOperatorOptions(t)

                # Add the clusters
                basename = conf['CLUSTER_DIR']+"/"
                add_clusters_with_colortable(basename,"_cluster_","cluster_colors",num_colors)

                # Add modes as labels
                visitUtils.add_labels(label_file,"geometrical_center")
                
                # save as image
                DrawPlots()
                
                visitUtils.save_window("tracking_",1)
                image_count=image_count+1;
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
        
        # clean up
        DeleteAllPlots();
        ClearWindow()
        
        if source_open:
            CloseDatabase(netcdf_file)
        CloseDatabase(label_file)
        visitUtils.close_pattern(basename+"*.vtr")
        visitUtils.close_pattern(basename+"*.vtk")
        return_code=call("rm -f *.vt*", shell=True)
        
        # periodically kill computing engine to
        # work around the memory leak fix
        image_count=image_count+1;
        if image_count % 100 == 0:
            CloseComputeEngine()

    # Close mapstuff
    close_topography()

    # create loops
    if conf['CREATE_SOURCE_MOVIE']:
        visitUtils.create_movie("source_","source.gif")
        visitUtils.create_movie("source_","source.m4v")

    if conf['CREATE_CLUSTERS_MOVIE']:
        visitUtils.create_movie("tracking_","tracking.gif")
        visitUtils.create_movie("tracking_","tracking.m4v")

    # clean up
    print "Cleaning up ..."
    return_code=call("mkdir images", shell=True)
    return_code=call("mv *.png images", shell=True)
    return_code=call("mkdir movies", shell=True)
    return_code=call("mv *.gif *.m4v movies", shell=True)
    return_code=call("rm -f *.vt* visitlog.py", shell=True)

    return


