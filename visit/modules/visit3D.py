#!/usr/bin/python
# Filename: visit3D.py

version = 'v1.4'

from visit import *
import sys
sys.path.append(".")
import glob
import os
import shutil
import time
import visitUtils
from subprocess import call

MAPS_HOME=os.environ['MEANIE3D_MAPDATA']
MAPSTUFF_FILE=MAPS_HOME+"/oase-mapdata.nc"
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

# Add 2D local topography, rivers and borders
# @param "local" or "national"
#
def add_mapstuff(extent):
    
    # Topography
    add_map_topography(extent);
    
    # Rivers & Boundaries
    add_map_borders(extent);
    add_map_rivers(extent);
    
    return

# Add 2D local topography
# @param "local" or "national"
#
def add_map_topography(extent):
    
    # open the file and add the plot
    OpenDatabase(MAPSTUFF_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
    # Topography
    AddPlot("Pseudocolor", extent+"_topo_3D")
    p = PseudocolorAttributes()
    p.pointSizePixels = 2
    p.invertColorTable=1
    p.colorTableName = "xray"
    p.scaling = p.Skew
    p.skewFactor=0.0001
    
    p.lightingFlag = 0
    p.legendFlag = 0;
    p.opacity = 1
    SetPlotOptions(p)
    
    return


# Add 2D borders
# @param "local" or "national"
#
def add_map_borders(extent):
    
    # open the file and add the plot
    OpenDatabase(MAPSTUFF_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
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
    
    return

# Add 2D borders
# @param "local" or "national"
#
def add_map_rivers(extent):
    
    # open the file and add the plot
    OpenDatabase(MAPSTUFF_FILE)
    
    if extent!="local" and extent !="national":
        print "ERROR:add_backdrop only accepts 'local' or 'national' as argument"
        return
    
    # Rivers & Boundaries
    
    AddPlot("Pseudocolor", "as_zonal/"+extent+"_rivers_3D")
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

# Closes databases connected with topo
def close_mapstuff():
    CloseDatabase(MAPSTUFF_FILE);
    return

#
# Sets default 3D view params for RADOLAN grid
# @param "local" or "national"
# @param perspective 1,2
# @param scale factor
#
def set_view_to_radolan(extend,perspective,scale_factor_z):

    v = GetView3D();
    
    if extend == "local":
    
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

    elif extend == "national":
    
        if perspective==1:
            v.viewNormal = (0.0244371, -0.668218, 0.743564)
            v.focus = (-73.9622, -4209.15, 7.31354)
            v.viewUp = (0.00033399, 0.743792, 0.668411)
            v.viewAngle = 30
            v.parallelScale = 636.44
            v.nearPlane = -1272.88
            v.farPlane = 1272.88
            v.imagePan = (0.00341995, 0.049739)
            v.imageZoom = 1.21
            v.perspective = 1
            v.eyeAngle = 2
            v.centerOfRotationSet = 0
            v.centerOfRotation = (0, 0, 0)
            v.axis3DScaleFlag = 0
            v.axis3DScales = (1, 1, 1)
            v.shear = (0, 0, 1)
        
        elif perspective==2:
            v.viewNormal = (-0.0310283, -0.976439, 0.213551)
            v.focus = (-73.9622, -4209.15, 7.31354)
            v.viewUp = (0.00253423, 0.213576, 0.976923)
            v.viewAngle = 30
            v.parallelScale = 636.44
            v.nearPlane = -1272.88
            v.farPlane = 1272.88
            v.imagePan = (0.0182949, -0.0442553)
            v.imageZoom = 1.4641
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
    return

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
# 'SCALE_FACTOR_Z' : scale factor for the Z axis.
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
def create_multiperspective(conf):

    DYLD_LIBRARY_PATH="/usr/local/lib"
    bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
    conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"

    # Set view and annotation attributes

    print "Setting annotation attributes:"
    set_annotations()
    
    if conf['RESUME'] == False:
        print "Removing results from previous runs"
        shutil.rmtree('images')
        shutil.rmtree('movies')
        return_code=call("rm -f *.vtk *.vtr *.png", shell=True)
    else:
        print "Removing intermediary files from previous runs"
        return_code=call("rm -f *.vtk *.vtr", shell=True)

    # Setting 3D view parameters
    print "Setting 3D view parameters"
    set_view_to_radolan(conf['GRID_EXTENT'],1,conf['SCALE_FACTOR_Z']);

    print "Creating colortables"
    num_colors = visitUtils.create_cluster_colortable("cluster_colors")

    if conf['WITH_TOPOGRAPHY']:
        visitUtils.create_topography_colortable()

    if conf['WITH_BACKGROUND_GRADIENT']:
        print "Setting background gradient"
        visitUtils.add_background_gradient();

    # Glob the netcdf directory
    print "Processing files in directory " + conf['NETCDF_DIR']
    netcdf_files = glob.glob(conf['NETCDF_DIR']+"/*.nc");

    # Keep track of number of images to allow
    # forced re-set in time to circumvent the
    # Visit memory leak
    image_count=0

    for netcdf_file in netcdf_files:
        
        # construct the cluster filename and find it
        # in the cluster directory
        
        netcdf_path,filename    = os.path.split(netcdf_file);
        basename                = os.path.splitext(filename)[0]
        cluster_file            = conf['CLUSTER_DIR']+"/"+basename+"-clusters.nc"
        label_file              = basename+"-clusters_centers.vtk"
        
        print "netcdf_file  = " + netcdf_file
        print "cluster_file = " + cluster_file
        
        # check if the files both exist
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
            
            OpenDatabase(netcdf_file);
            source_open = True

            if conf['WITH_TOPOGRAPHY']:
                # add 3D topograpy
                add_mapstuff("local")

            if conf['CREATE_SOURCE_MOVIE']:
                
                if conf['WITH_TOPOGRAPHY']:
                    print "-- Adding topography data --"
                    add_topography("national_topo_3D")
                
                if conf['WITH_RIVERS_AND_BOUNDARIES']:
                    print "-- Adding map data --"
                    add_map_rivers("national")
                    add_map_borders("national")
                
                if conf['WITH_DATETIME']:
                    print "-- Adding timestamp --"
                    visitUtils.add_datetime(netcdf_file)
                
                print "-- Plotting source data --"
                start_time = time.time()
                
                # Add source data and threshold it

                variables = conf['VARIABLES']

                for i in range(len(variables)):

                    add_pseudocolor(netcdf_file,
                                            conf['VARIABLES'][i],
                                            conf['COLORTABLES'][i],
                                            conf['OPACITY'][i],1)
                    p = PseudocolorAttributes()
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
                f1 = "p1_tracking_" + number_postfix
                f2 = "p2_tracking_" + number_postfix
                                 
                if os.path.exists(f1) and os.path.exists(f2):
                    print "Skipping existing files " + f1 + "," + f2
                    skip = True
                elif (os.path.exists(f1) and not os.path.exists(f2)):
                    os.remove(f1)
                elif (os.path.exists(f2) and not os.path.exists(f1)):
                    os.remove(f2)
                                 
            if skip == False:
 
                start_time = time.time()
                print "-- Converting clusters to .vtr --"
                
                # build the clustering command
                command=conversion_bin+" -f "+cluster_file+" "+conf['CONVERSION_PARAMS']
                print command
                return_code = call( command, shell=True)
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
                
                print "-- Rendering cluster scene --"
                start_time = time.time()
                
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
                
                DrawPlots()
                
                # save as image
                visitUtils.save_window("p1_tracking_",1)
                
                # change perspective
                set_view_to_radolan(conf['GRID_EXTENT'],2,conf['SCALE_FACTOR_Z']);
                
                visitUtils.save_window("p2_tracking_",1)
                
                # change perspective back
                set_view_to_radolan(conf['GRID_EXTENT'],1,conf['SCALE_FACTOR_Z']);
                
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

    # close mapstuff

    close_mapstuff()

    # create loops

    if conf['CREATE_SOURCE_MOVIE']:
        visitUtils.create_movie("source_","source.gif")
        visitUtils.create_movie("source_","source.m4v")

    if  conf['CREATE_CLUSTERS_MOVIE']:
        visitUtils.create_movie("p1_tracking_","p1_tracking.gif")
        visitUtils.create_movie("p1_tracking_","p1_tracking.m4v")
        visitUtils.create_movie("p2_tracking_","p2_tracking.gif")
        visitUtils.create_movie("p2_tracking_","p2_tracking.m4v")

    # clean up

    print "Cleaning up ..."
    return_code=call("mkdir images", shell=True)
    return_code=call("mv *.png images", shell=True)
    return_code=call("mkdir movies", shell=True)
    return_code=call("mv *.gif *.m4v movies", shell=True)
    return_code=call("rm -f *.vt* visitlog.py", shell=True)

    return

