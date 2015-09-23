#!/usr/bin/python

# ------------------------------------------------------------------------------
# Filename: visit3D.py
#
# This module bundles python routines for handling Visit3D plotting
#
# @author Juergen Simon (juergen_simon@mac.com)
# ------------------------------------------------------------------------------

from visit import *
import sys
sys.path.append(".")

import glob
import os
import time
import meanie3D_visit_utils
from pprint import pprint
from subprocess import call

# ------------------------------------------------------------------------------
# The following variables are important for referencing the mapstuff and
# topography data
# ------------------------------------------------------------------------------
MAPS_HOME=os.environ['MEANIE3D_MAPDATA']
MAPSTUFF_FILE=MAPS_HOME+"/oase-mapdata.nc"
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

# ------------------------------------------------------------------------------
# Adds clusters with names "*_infix_*.vtk" to the
# current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
# ------------------------------------------------------------------------------
def add_clusters_with_colortable(basename,infix,color_table_name,color_count,conf):
    
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
            
        cluster_opacity = 1.0
        if conf["CLUSTER_OPACITY"]:
            cluster_opacity = conf["CLUSTER_OPACITY"]
            
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
        # cp.opacityType = cp.Constant
        cp.opacity=cluster_opacity
        SetPlotOptions(cp)
        
    return

# ------------------------------------------------------------------------------
# Adds boundaries with to the current visit window.
# @param basename first part of the search pattern
#        used to find the clusters
# @param infix (e.g. "_untracked_clusters_"
# @param color tables
# @param number of colors in the color table
# ------------------------------------------------------------------------------
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
    p.pointSizePixels = 2
    p.colorTableName = color_table_name
    p.lightingFlag = 0
    p.legendFlag = legendFlag;
    p.opacity = opacity
    SetPlotOptions(p)
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

# ------------------------------------------------------------------------------
# Add 2D borders
# @param "local" or "national"
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Add 2D borders
# @param "local" or "national"
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Closes databases connected with topo
# ------------------------------------------------------------------------------


def add_displacement_vectors(displacements_file):
    # open displacements vector file and add plot
    OpenDatabase(displacements_file)
    AddPlot("Vector","displacement")
    # Set plot's attributes
    p=VectorAttributes();
    p.useStride = 1
    p.stride = 1
    p.scale = 0.25
    p.scaleByMagnitude = 1
    p.autoScale = 0
    p.headSize = 0.25
    p.headOn = 1
    p.useLegend = 0
    SetPlotOptions(p)
    return
    
def close_mapstuff():
    CloseDatabase(MAPSTUFF_FILE);
    return

# ------------------------------------------------------------------------------
# Sets the view up with the given perspective object
# @param perspective object
# @param stretching factor for z axis
# ------------------------------------------------------------------------------

def set_perspective(perspective,scale_factor_z):
    # get a view config obect
    v = GetView3D();
    
    # set the values from the json perspective
    # dictionary
    
    v.viewNormal = tuple(perspective["viewNormal"])
    v.focus = tuple(perspective["focus"])
    v.viewUp = tuple(perspective["viewUp"])
    v.viewAngle = perspective["viewAngle"]
    v.parallelScale = perspective["parallelScale"]
    v.nearPlane = perspective["nearPlane"]
    v.farPlane = perspective["farPlane"]
    v.imagePan = tuple(perspective["imagePan"])
    v.imageZoom = perspective["imageZoom"]
    v.perspective = perspective["perspective"]
    v.eyeAngle = perspective["eyeAngle"]
    v.centerOfRotationSet = perspective["centerOfRotationSet"]
    v.centerOfRotation = tuple(perspective["centerOfRotation"])
    v.axis3DScaleFlag = perspective["axis3DScaleFlag"]
    v.axis3DScales = tuple(perspective["axis3DScales"])
    v.shear = tuple(perspective["shear"])

    # scale z axis
    
    v.axis3DScaleFlag = 1
    v.axis3DScales = (1, 1, scale_factor_z)
    
    # apply
    
    SetView3D(v)
    return

# ------------------------------------------------------------------------------
# Sets default 3D view params for RADOLAN grid
# @param "local" or "national"
# @param stretching factor for z axis
# ------------------------------------------------------------------------------
def set_view_to_radolan(extend,scale_factor_z):

    v = GetView3D();
    
    if extend == "local":
    
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
        
    elif extend == "national":
    
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
        
    if scale_factor_z != 1.0:
        v.axis3DScaleFlag = 1
        v.axis3DScales = (1, 1, scale_factor_z)

    SetView3D(v);
    return

# ------------------------------------------------------------------------------
# Sets up standard values for axis etc
# ------------------------------------------------------------------------------
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
# Checks if the file with the given number exists. Takes
# perspectives into account. If any perspective exists,
# all images from this number are deleted
# @param configuration
# @param basename ('source_','tracking_' etc.)
# @param number number of file
# ------------------------------------------------------------------------------
def delete_images(conf,basename,image_count):
    number_postfix = str(image_count).rjust(4,'0') + ".png";
    result = False
    if 'PERSPECTIVES' in conf.keys():
        perspective_nr = 1
        for perspective in conf['PERSPECTIVES']:
            fn = "p"+str(perspective_nr)+"_"+basename+"_"+number_postfix
            if (os.path.exists(fn)):
                os.remove(fn)
            perspective_nr = perspective_nr + 1
    else:
        fn = basename+"_"+number_postfix
        if (os.path.exists(fn)):
            os.remove(fn)

# ------------------------------------------------------------------------------
# Checks if the file with the given number exists. Takes
# perspectives into account
# @param configuration
# @param basename ('source_','tracking_' etc.)
# @param number number of file
# @return "all","none" or "partial"
# ------------------------------------------------------------------------------
def images_exist(conf,basename,image_count):
    number_postfix = str(image_count).rjust(4,'0') + ".png";
    result = False
    if 'PERSPECTIVES' in conf.keys():
        num_perspectives = len(conf['PERSPECTIVES'])
        num_found = 0
        perspective_nr = 1
        for perspective in conf['PERSPECTIVES']:
            fn = "p"+str(perspective_nr)+"_"+basename+"_"+number_postfix
            if (os.path.exists(fn)):
                num_found = num_found + 1
            perspective_nr = perspective_nr + 1
        if num_found == 0: 
            return "none"
        elif num_found == num_perspectives:
            return "all"
        else:
            return "partial"
    else:
        fn = basename+"_"+number_postfix
        if os.path.exists(fn):
            return "all"
        else:
            return "none"

# ------------------------------------------------------------------------------
# Generic routine for visualizing 3D clusters in two perspectives
#
# The following configuration options exist:
#
# 'source_directory' : directory with the source data files
# 'cluster_directory' : directory with the cluster results
# 'meanie3d_home' : home directory of meanie3D (for the mapstuff file and modules)
# 'resume' : if true, the existing image files are not wiped and work is
#            picked up where it left off. Otherwise all existing images
#            are deleted and things are started from scratch
# 'with_background_gradient' : add a gray background gradient to the canvas?
# 'with_topography' : use the topography data from the mapstuff file?
# 'with_rivers_and_boundaries' : add rivers and boundaries?
# 'with_source_backround' : re-add the source data when plotting clusters?
# 'with_datetime' : add a date/time label?
# 'create_source_movie' : create a movie from the source images?
# 'create_clusters_movie' : create a movie from the cluster images?
# 'SCALE_FACTOR_Z' : scale factor for the Z axis.
# 'grid_extent' : "national" or "local"
# 'conversion_params' : parameters for meanie3D-cfm2vtk
# 'variables' : list of variables for the source data
# 'lower_tresholds' : bottom cutoff for each variable,
# 'upper_tresholds' : top cutoff for each variable,
# 'var_min' : lowest value on legend
# 'var_max' : highest value on legend
# 'colortables' : colortable to use for each variable
# 'colortables_invert_flags' : flag indicating inversion of colortable
#                              for each variable
# 'PERSPECTIVES' : array with perspective objects
# 'opacity' : opacity to use for each variable
# ------------------------------------------------------------------------------
def visualization(conf):
    
    print "-------------------------------------------------"
    print "visit3D.visualization started with configuration:"
    print "-------------------------------------------------"
    pprint(conf)
    print "-------------------------------------------------"

    bin_prefix = "export DYLD_LIBRARY_PATH="+meanie3D_visit_utils.get_dyld_library_path()+";"
    conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"

    # Set view and annotation attributes

    print "Setting annotation attributes:"
    set_annotations()
    
    if conf['resume'] == False:
        print "Removing results from previous runs"
        return_code=call("rm -rf images movies *.vtk *.vtr *.png", shell=True)
    else:
        print "Removing intermediary files from previous runs"
        return_code=call("rm -f *.vtk *.vtr", shell=True)

    if conf['create_clusters_movie']:
        print "Creating colortables"
        num_colors = meanie3D_visit_utils.create_cluster_colortable("cluster_colors")

    if conf['with_topography']:
        meanie3D_visit_utils.create_topography_colortable()

    if conf['with_background_gradient']:
        print "Setting background gradient"
        meanie3D_visit_utils.add_background_gradient();

    scaleFactorZ = 1.0
    if "SCALE_FACTOR_Z" in conf.keys():
        scaleFactorZ = conf['SCALE_FACTOR_Z']

    # Glob the netcdf directory
    print "Processing files in directory " + conf['source_directory']
    netcdf_files = glob.glob(conf['source_directory']+"/*.nc");

    # Keep track of number of images to allow
    # forced re-set in time to circumvent the
    # Visit memory leak
    image_count=0

    for netcdf_file in netcdf_files:
        
        # construct the cluster filename and find it
        # in the cluster directory
        
        netcdf_path,filename    = os.path.split(netcdf_file);
        basename                = os.path.splitext(filename)[0]
        cluster_file            = conf['cluster_directory']+"/"+basename+"-clusters.nc"
        label_file              = basename+"-clusters_centers.vtk"
        displacements_file      = basename+"-clusters_displacements.vtk"
        
        print "netcdf_file  = " + netcdf_file
        print "cluster_file = " + cluster_file
        
        # check if the files both exist
        if not os.path.exists(cluster_file):
            print "Cluster file does not exist. Skipping."
            continue

        # predict the filenames for checking on resume
        number_postfix = str(image_count).rjust(4,'0')

        source_open = False
        skip_source = False
        
        if conf['resume'] == True:
            exists = images_exist(conf,"source",image_count)
            if exists == "all":
                print "Source visualization "+number_postfix+" exists. Skipping."
                skip_source = True
            elif exists == "partial":
                print "Deleting partial visualization " + number_postfix
                delete_images(conf,"source",image_count)

        if skip_source == False:
            
            OpenDatabase(netcdf_file);
            source_open = True

            if conf['create_source_movie']:
                
                if conf['with_topography']:
                    print "-- Adding topography data --"
                    add_map_topography(conf['grid_extent'])
                
                if conf['with_rivers_and_boundaries']:
                    print "-- Adding map data --"
                    add_map_rivers(conf['grid_extent'])
                    add_map_borders(conf['grid_extent'])
                
                if conf['with_datetime']:
                    print "-- Adding timestamp --"
                    meanie3D_visit_utils.add_datetime(netcdf_file)
                
                print "-- Plotting source data --"
                start_time = time.time()
                
                # Add source data and threshold it

                variables = conf['variables']

                for i in range(len(variables)):

                    add_pseudocolor(netcdf_file,
                                            conf['variables'][i],
                                            conf['colortables'][i],
                                            conf['opacity'][i],1)
                    p = PseudocolorAttributes()
                    p.minFlag,p.maxFlag=1,1
                    p.min,p.max=conf['var_min'][i],conf['var_max'][i]
                    p.invertColorTable = conf['colortables_invert_flags'][i]
                    SetPlotOptions(p)
                
                    AddOperator("Threshold")
                    t = ThresholdAttributes();
                    t.lowerBounds=(conf['lower_tresholds'][i])
                    t.upperBounds=(conf['upper_tresholds'][i])
                    SetOperatorOptions(t)

                DrawPlots();
                
                if 'PERSPECTIVES' in conf.keys():
                    perspective_nr = 1
                    for perspective in conf['PERSPECTIVES']:
                        set_perspective(perspective,scaleFactorZ)
                        filename = "p" + str(perspective_nr) + "_source_"
                        meanie3D_visit_utils.save_window(filename,1)
                        perspective_nr = perspective_nr + 1
                else:
                    set_view_to_radolan(conf['grid_extent'],conf['SCALE_FACTOR_Z'])
                    meanie3D_visit_utils.save_window("source_",1)
                
                DeleteAllPlots()
                ClearWindow()
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
        
        if conf['create_clusters_movie']:
                                 
            skip = False
            
            if conf['resume'] == True:

                exists = images_exist(conf,"tracking",image_count)
                if exists == "all":
                    print "Cluster visualization "+number_postfix+" exists. Skipping."
                    skip = True
                elif exists == "partial":
                    print "Deleting partial cluster visualization " + number_postfix
                    delete_images(conf,"tracking",image_count)
                                 
            if skip == False:
 
                start_time = time.time()
                print "-- Converting clusters to .vtr --"
                
                # build the clustering command
                command=conversion_bin+" -f "+cluster_file+" "+conf['conversion_params']
                if conf['WITH_DISPLACEMENT_VECTORS']:
                    command = command + " --write-displacement-vectors"
                print command
                return_code = call( command, shell=True)
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)
                
                print "-- Rendering cluster scene --"
                start_time = time.time()
                
                if conf['with_topography']:
                    print "-- Adding topography data --"
                    add_map_topography(conf['grid_extent'])
                
                if conf['with_rivers_and_boundaries']:
                    print "-- Adding map data --"
                    add_map_rivers(conf['grid_extent'])
                    add_map_borders(conf['grid_extent'])
                
                if conf['with_datetime']:
                    print "-- Adding timestamp --"
                    meanie3D_visit_utils.add_datetime(netcdf_file)
                
                if conf['with_source_backround']:
                    
                    if not source_open:
                        OpenDatabase(netcdf_file);
                        source_open = True

                    for i in range(len(variables)):
                    
                        add_pseudocolor(netcdf_file,
                                                conf['variables'][i],
                                                conf['colortables'][i],
                                                conf['opacity'][i],1)
                        p = PseudocolorAttributes()
                        p.invertColorTable = conf['colortables_invert_flags'][i]
                        p.minFlag,p.maxFlag=1,1
                        p.min,p.max=conf['var_min'][i],conf['var_max'][i]
                        SetPlotOptions(p)
                        
                        AddOperator("Threshold")
                        t = ThresholdAttributes();
                        t.lowerBounds=(conf['lower_tresholds'][i])
                        t.upperBounds=(conf['upper_tresholds'][i])
                        SetOperatorOptions(t)

                
                # Add the clusters
                basename = conf['cluster_directory']+"/"
                add_clusters_with_colortable(basename, "_cluster_", "cluster_colors", num_colors, conf)
                
                # Add modes as labels
                meanie3D_visit_utils.add_labels(label_file,"geometrical_center")
                
                # Add displacement vectors
                if conf['WITH_DISPLACEMENT_VECTORS']:
                    add_displacement_vectors(displacements_file)
                
                DrawPlots()
                
                if 'PERSPECTIVES' in conf.keys():
                    perspective_nr = 1
                    for perspective in conf['PERSPECTIVES']:
                        set_perspective(perspective,scaleFactorZ)
                        filename = "p" + str(perspective_nr) + "_tracking_"
                        meanie3D_visit_utils.save_window(filename,1)
                        perspective_nr = perspective_nr + 1
                else:
                    set_view_to_radolan(conf['grid_extent'],conf['SCALE_FACTOR_Z'])
                    meanie3D_visit_utils.save_window("tracking_",1)
                
                # change perspective back
                set_view_to_radolan(conf['grid_extent'],conf['SCALE_FACTOR_Z']);
                
                print "    done. (%.2f seconds)" % (time.time()-start_time)

        # clean up
        
        DeleteAllPlots();
        ClearWindow()
        if source_open:
            CloseDatabase(netcdf_file)
            CloseDatabase(label_file)
        meanie3D_visit_utils.close_pattern(basename+"*.vtr")
        meanie3D_visit_utils.close_pattern(basename+"*.vtk")
        return_code=call("rm -f *.vt*", shell=True)
        
        # periodically kill computing engine to
        # work around the memory leak fix
        image_count=image_count+1;
        if image_count % 100 == 0:
            CloseComputeEngine()

    # close mapstuff

    close_mapstuff()

    # create loops

    if 'PERSPECTIVES' in conf.keys():

        perspective_nr = 1

        for perspective in conf['PERSPECTIVES']:

            if conf['create_source_movie']:
                movie_fn = "p" + str(perspective_nr) + "_source"
                image_fn = movie_fn + "_"
                meanie3D_visit_utils.create_movie(image_fn,movie_fn+".gif")
                meanie3D_visit_utils.create_movie(image_fn,movie_fn+".m4v")

            if  conf['create_clusters_movie']:
                movie_fn = "p" + str(perspective_nr) + "_tracking"
                image_fn = movie_fn + "_"
                meanie3D_visit_utils.create_movie(image_fn,movie_fn+".gif")
                meanie3D_visit_utils.create_movie(image_fn,movie_fn+".m4v")

            perspective_nr = perspective_nr + 1
    else:

        if conf['create_source_movie']:
            meanie3D_visit_utils.create_movie("source_","source.gif")
            meanie3D_visit_utils.create_movie("source_","source.m4v")

        if  conf['create_clusters_movie']:
            meanie3D_visit_utils.create_movie("tracking_","tracking.gif")
            meanie3D_visit_utils.create_movie("tracking_","tracking.m4v")

    # clean up

    print "Cleaning up ..."
    return_code=call("mkdir images", shell=True)
    return_code=call("mv *.png images", shell=True)
    return_code=call("mkdir movies", shell=True)
    return_code=call("mv *.gif *.m4v movies", shell=True)
    return_code=call("rm -f *.vt* visitlog.py", shell=True)

    return

