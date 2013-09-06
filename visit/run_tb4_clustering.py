#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
DYLD_LIBRARY_PATH = "DL_PATH"
NETCDF_DIR        = "SOURCE_DIR"
VAR_NAME          = "VAR_P"

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

CLUSTERING_PARAMS  = " -d lat,lon --vtk-dimensions lon,lat --ranges=1,1,300"
CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " --write-variables-as-vtk "+VAR_NAME
CLUSTERING_PARAMS += " --upper-thresholds "+VAR_NAME+"=255"
CLUSTERING_PARAMS += " --weight-function inverse"
CLUSTERING_PARAMS += " --write-clusters-as-vtk"
CLUSTERING_PARAMS += " --write-cluster-weight-response"
CLUSTERING_PARAMS += " --verbosity 1"

#TRACKING_PARAMS = "--verbosity 1 --write-vtk"
#TRACKING_PARAMS += " --vtk-dimensions x,y"
#TRACKING_PARAMS += " -t "+VAR_NAME
#TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

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

#Cluster color tables
#col_tables = ["Purples","Blues","Oranges","Greens","Reds"]

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

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)
last_cluster_file=""

# Set up viewport
v = GetView2D()
v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
v.windowCoords = (3, 19.8, 45, 55.8)
SetView2D(v)

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:
    
    basename = os.path.basename(netcdf_file)
    stripped_name=os.path.splitext(basename)[0];
    cluster_file=stripped_name+"-clusters.nc"
    vtk_file=stripped_name + "_" + VAR_NAME + ".vtk"

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

    # the first file in the specific series caused some trouble
    # with the viewport. Comment this out as the data becomes
    # more consistent

    visit2D.add_pseudocolor( netcdf_file, VAR_NAME, "hot_desaturated",1 )
    cp=PseudocolorAttributes();
    cp.invertColorTable=1
    SetPlotOptions(cp)

    DrawPlots()
    RecenterView()

    # Calling ToggleMaintainViewMode helps
    # keeping the window from 'jittering'
    #if toggled_maintain != True :
        #ToggleMaintainViewMode()
        #toggled_maintain=True
    
    #visitUtils.save_window(stripped_name+"_source_",1)
    
    print "    done. (%.2f seconds)" % (time.time()-start_time)
    print "-- Rendering clusters --"
    start_time = time.time()

    # plot the clusters
    #cluster_basename = "Release/*_cluster_*.vtk"

    cluster_basename = os.path.splitext(basename)[0] + "*-clusters_weight*.vtk"
    list = glob.glob( cluster_basename )

    count = 0;
    for fname in list:
    
        # add plot
        OpenDatabase(fname);
        AddPlot("Pseudocolor", "weight")
        
        # set plot attributes accordingly
        cp=PseudocolorAttributes();
        cp.pointSizePixels=10
        cp.opacity=0.25
        cp.minFlag,cp.maxFlag=1,1
        cp.min,cp.max=0.0,0.5
        cp.legendFlag=0
        cp.colorTableName = "contoured";
        if count%2==0:
            cp.invertColorTable=1
        else:
            cp.invertColorTable=0
        SetPlotOptions(cp)
        count = count+1

    DrawPlots()
    RecenterView()
    visitUtils.save_window(stripped_name+"_clusters_",1)

    #
    # clean up
    #

    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(netcdf_file)
    visitUtils.close_pattern(os.path.splitext(basename)[0]+"*-clusters_weight*.vtk")
    return_code=call("rm -f *cluster*.vtk", shell=True)
    return_code=call("rm -f *.vtk", shell=True)

    print "    done. (%.2f seconds)" % (time.time()-start_time)

    # keep track
    last_cluster_file=cluster_file

    # don't forget to increment run counter
    run_count = (run_count + 1)

print "Done. Closing Visit."
exit()

