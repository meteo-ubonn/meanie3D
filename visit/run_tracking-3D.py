#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
DYLD_LIBRARY_PATH = "DL_PATH"
NETCDF_DIR        = "SOURCE_DIR"
PARAM_T           = "SCALE"

# Appending the module path is crucial

sys.path.append(MEANIE3D_HOME+"/visit/modules")

import glob
import os
import sys
import time
import visit3D
import visitUtils
from subprocess import call

print [key for key in locals().keys()
       if isinstance(locals()[key], type(sys)) and not key.startswith('__')]

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# RADOLAN

VAR_NAME="zh"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds zh=30 -m 10"

CLUSTERING_PARAMS =  "-d z,y,x --vtk-dimensions x,y,z"
CLUSTERING_PARAMS += " --verbosity 1"
CLUSTERING_PARAMS += " --write-clusters-as-vtk"
CLUSTERING_PARAMS += " --write-variables-as-vtk="+VAR_NAME
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS


TRACKING_PARAMS = "--verbosity 1 --write-vtk"
TRACKING_PARAMS += " --vtk-dimensions x,y,z"
TRACKING_PARAMS += " -t "+VAR_NAME
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME
print "DYLD_LIBRARY_PATH="+DYLD_LIBRARY_PATH

toggled_maintain=False

# delete previous results
if last_completed_run_count == 0:
    print "Cleaning up *.vtk *.nc *.png *.log"
    return_code=call("rm -f *.nc *.vtk *.log *.png", shell=True)

# binaries
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"

#detection_bin = bin_prefix + "M3D_HOME/Release/" + "meanie3D-detect"
#tracking_bin  = bin_prefix + "M3D_HOME/Release/" + "meanie3D-track"
#trackplot_bin = bin_prefix + "M3D_HOME/Release/" + "meanie3D-trackplot"

detection_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-detect"
tracking_bin  = bin_prefix + "/usr/local/bin/" + "meanie3D-track"
trackplot_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-trackplot"

print "Detection Command:"
print detection_bin
print "Tracking Command:"
print tracking_bin
print "Tracking Evaluation Command:"
print trackplot_bin

# Cluster color tables
col_tables = ["Purples","Blues","Oranges","Greens","Reds","Paired"]

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes
a = GetAnnotationAttributes()
a.axes2D.visible=1
a.axes2D.autoSetScaling=0
a.userInfoFlag=0
a.timeInfoFlag=0
a.legendInfoFlag=0
a.databaseInfoFlag=1
SetAnnotationAttributes(a)

# Modify view parameters
#v = GetView2D()
#v.windowCoords = (-418.462, 292.538, -4446.64, -3759.64)
#v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
#SetView2D(v)

v = GetView3D();
v.viewNormal = (0.656802,-0.498223,0.566025)
v.focus = (-239.212,-4222.9,7.375)
v.viewUp = (-0.457525,0.333371,0.824339)
v.viewAngle = 30
v.parallelScale = 173.528
v.nearPlane = -347.056
v.farPlane = 347.056
v.imagePan = (0, 0)
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

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)
last_cluster_file=""

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:
    
    basename = os.path.basename(netcdf_file)
    cluster_file=os.path.splitext(basename)[0]+"-clusters.nc"
    vtk_file=os.path.splitext(basename)[0] + "_" + VAR_NAME + ".vtk"
    
    # if there is a resume counter, keep skipping
    # until the count is right
    if (last_completed_run_count > 0) and (run_count <= last_completed_run_count):
        last_cluster_file=cluster_file
        run_count = run_count + 1
        continue
    
    print "----------------------------------------------------------"
    print "Processing " + netcdf_file
    print "----------------------------------------------------------"
    
    print "-- Clustering --"
    start_time = time.time()
    
    #
    # Cluster
    #
    
    # build the clustering command
    command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + CLUSTERING_PARAMS

    # use previous result to enhance current
    if run_count > 0:
        command = command + " -p " + last_cluster_file

    command = command + " > clustering_" + str(run_count)+".log"

    # execute
    print command
    return_code = call( command, shell=True)
    
    print "    done. (%.2f seconds)" % (time.time()-start_time)
    print "-- Rendering source data --"
    start_time = time.time()
    
    #
    # Plot the source data in color
    #
    
    visit3D.add_pseudocolor( vtk_file, VAR_NAME, "hot_desaturated",0.25 )
    DrawPlots()
    
    # Calling ToggleMaintainViewMode helps
    # keeping the window from 'jittering'
    if toggled_maintain != True :
        ToggleMaintainViewMode()
        toggled_maintain=True
    
    visitUtils.save_window("source_",1)
    
    # clean up
    DeleteAllPlots();
    
    print "    done. (%.2f seconds)" % (time.time()-start_time)
    print "-- Rendering untracked clusters --"
    start_time = time.time()
    
    #
    # Plot untracked clusters
    #
    
    # Re-add the source with "xray"
    # visit3D.add_pseudocolor(vtk_file,VAR_NAME,"xray",0)
    
    # Add the clusters
    visit3D.add_clusters(basename,"_cluster_",col_tables)
    
    # Add modes as labels
    label_file=os.path.splitext(basename)[0]+"-clusters_centers.vtk"
    visitUtils.add_labels(label_file,"geometrical_center")
    
    # Get it all processed and stowed away
    DrawPlots()
    visitUtils.save_window("untracked_",1)
    
    #
    # clean up
    #
    
    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(vtk_file)
    CloseDatabase(label_file)
    visitUtils.close_pattern(basename+"*_cluster_*.vtk")
    return_code=call("rm -f *cluster*_*.vtk", shell=True)
    
    print "    done. (%.2f seconds)" % (time.time()-start_time)
    
    #
    # Tracking
    #
    
    # if we have a previous scan, run the tracking command
    
    if run_count > 0:
        
        print "-- Tracking --"
        start_time = time.time()
        
        command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + TRACKING_PARAMS
        command = command + " > tracking_" + str(run_count)+".log"
        
        # execute
        return_code = call( command, shell=True)
        
        print "    done. (%.2f seconds)" % (time.time()-start_time)
    
    print "-- Rendering tracked clusters --"
    start_time = time.time()
    
    #
    # Plot tracked clusters
    #
    
    # Re-add the source with "xray"
    visit3D.add_pseudocolor(vtk_file,VAR_NAME,"xray",0)
    
    if (run_count > 0):
        
        # Add the clusters
        visit3D.add_clusters(basename,"_cluster_",col_tables)
        
        # Add modes as labels
        label_file=os.path.splitext(basename)[0]+"-clusters_centers.vtk"
        visitUtils.add_labels(label_file,"geometrical_center")
    
    # Get it all processed and stowed away
    DrawPlots()
    visitUtils.save_window("tracked_",1)

    #
    # clean up
    #
    
    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(vtk_file)
    
    if run_count > 0:
        CloseDatabase(label_file)
        visitUtils.close_pattern(basename+"*_cluster_*.vtk")
        return_code=call("rm -f *cluster_*.vtk", shell=True)
    
    print "    done. (%.2f seconds)" % (time.time()-start_time)
    
    # keep track
    last_cluster_file=cluster_file
    
    # don't forget to increment run counter
    run_count = (run_count + 1)


print "Done. Closing Visit."
exit()
