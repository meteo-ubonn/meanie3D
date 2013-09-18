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
import visit2D
import visitUtils
from subprocess import call

print [key for key in locals().keys()
       if isinstance(locals()[key], type(sys)) and not key.startswith('__')]

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# RADOLAN

VAR_NAME="RX"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds RX=35 -m 5"

CLUSTERING_PARAMS =  "-d y,x --vtk-dimensions x,y"
CLUSTERING_PARAMS += " --verbosity 1"
CLUSTERING_PARAMS += " --write-variables-as-vtk="+VAR_NAME
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " --write-clusters-as-vtk"

CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS


TRACKING_PARAMS = "--verbosity 1 --write-vtk"
TRACKING_PARAMS += " --vtk-dimensions x,y"
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
col_tables = ["Purples","Blues","Oranges","Greens","Reds"]

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
v = GetView2D()
v.windowCoords = (-418.462, 292.538, -4446.64, -3759.64)
v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
SetView2D(v)

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)
last_cluster_file=""

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:
    
    basename = os.path.basename(netcdf_file)
    cluster_file=os.path.splitext(basename)[0]+"-clusters.nc"
    #vtk_file=os.path.splitext(basename)[0] + "_" + VAR_NAME + ".vtk"
    vtk_file=os.path.splitext(basename)[0] + ".vtr"

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
    
    visit2D.add_pseudocolor( vtk_file, VAR_NAME, "hot_desaturated",1 )
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
    visit2D.add_pseudocolor(vtk_file,VAR_NAME,"xray",0)
    
    # Add the clusters
    visit2D.add_clusters(basename,"_cluster_",col_tables,0.5)
    
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
    visit2D.add_pseudocolor(vtk_file,VAR_NAME,"xray",0)
    p = PseudocolorAttributes()
    p.opacity=0.25
    SetPlotOptions(p)


    if (run_count > 0):

        # Add the clusters
        visit2D.add_clusters(basename,"_cluster_",col_tables,0.5)

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
        return_code=call("rm -f *.vtk", shell=True)

    print "    done. (%.2f seconds)" % (time.time()-start_time)

    # keep track
    last_cluster_file=cluster_file

    # don't forget to increment run counter
    run_count = (run_count + 1)


print "Done. Closing Visit."
exit()

