#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
DYLD_LIBRARY_PATH = "DL_PATH"
NETCDF_DIR        = "SOURCE_DIR"

# Appending the module path is crucial 
sys.path.append(MEANIE3D_HOME+"/visit/modules")

import glob
import sys
import os
import visit2D
from subprocess import call

# OASE Komposit
#CLUSTERING_PARAMS = "-d z,y,x -v zh -w zh -r 5,10,10,100 --drf-threshold 0.75 -s 16 -t 10 --write-variables-as-vtk=zh --vtk-dimensions=x,y,z" 
#TRACKING_PARAMS   = "-t zh --verbosity 3 --write-vtk --vtk-dimensions=x,y,z --wr=1.0 --ws=1.0 --wt=0.0"

# RADOLAN
CLUSTERING_PARAMS = "--write-variables-as-vtk=reflectivity -v reflectivity -w reflectivity -d x,y -r 5,5,200 --drf-threshold 0.5 -s 64 -t 20 -m 10 --verbosity 1 "
TRACKING_PARAMS   = "-t reflectivity --verbosity 1 --write-vtk --wr=1.0 --ws=0.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME
print "DYLD_LIBRARY_PATH="+DYLD_LIBRARY_PATH

toggled_maintain=False

# delete any previous results
# results are stored in the work directory

# delete previous results
return_code=call("rm *-clusters.nc", shell=True)
return_code=call("rm *_clusters_*.vtk", shell=True)

# binaries
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
detection_bin = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-detect"
tracking_bin  = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-track"

print "Detection Command:"
print detection_bin
print "Tracking Command:"
print tracking_bin

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
v = GetView3D()
v.focus=(-238.5,-4222.5,50.0)
SetView3D(v)

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)
last_cluster_file=""

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:

    print "----------------------------------------------------------"
    print "Processing " + netcdf_file
    print "----------------------------------------------------------"

    basename = os.path.basename(netcdf_file)
    cluster_file=os.path.splitext(basename)[0]+"-clusters.nc"
    vtk_file=os.path.splitext(basename)[0]+".vtk"

    print "-- Rendering source data --"

    #
    # Plot the source data in color
    #
    
    visit2D.add_pseudocolor_2D(vtk_file,"reflectivity","hot_desaturated")
    DrawPlots()
    
    # Calling ToggleMaintainViewMode helps
    # keeping the window from 'jittering'
    if toggled_maintain != True :
        ToggleMaintainViewMode()
        toggled_maintain=True
    
    visit2D.save_window("source_",1)
    
    # clean up
    DeleteAllPlots();

    print "-- Clustering --"

    #
    # Cluster
    #
    
    # build the clustering command
    command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + CLUSTERING_PARAMS
    command = command + " --write-clusters-as-vtk"
    command = command + " > clustering_" + str(run_count)+".log"

    # execute
    print command
    return_code = call( command, shell=True)
    
    print "-- Rendering untracked clusters --"

    #
    # Plot untracked clusters
    #
    
    # Re-add the source with "xray"
    visit2D.add_pseudocolor_2D(vtk_file,"reflectivity","xray")
    
    # Add the clusters
    visit2D.add_clusters_2D(basename,"_cluster_",col_tables)
    
    # Add modes as labels
    modes_file=os.path.splitext(basename)[0]+"-clusters_modes.vtk"
    visit2D.add_modes(modes_file)

    # Get it all processed and stowed away
    DrawPlots()
    visit2D.save_window("untracked_",1)

    #
    # clean up
    #

    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(vtk_file)
    CloseDatabase(modes_file)
    visit2D.close_clusters_2D(basename,"_cluster_")
    return_code=call("rm *cluster_*.vtk", shell=True)

    print "-- Tracking --"

    #
    # Tracking
    #

    # if we have a previous scan, run the tracking command

    if run_count > 0:
        command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + TRACKING_PARAMS
        command = command + " > tracking_" + str(run_count)+".log"
        
        # execute
        print command
        return_code = call( command, shell=True)

    # keep track
    last_cluster_file=cluster_file

    print "-- Rendering tracked clusters --"

    #
    # Plot tracked clusters
    #

    # Re-add the source with "xray"
    visit2D.add_pseudocolor_2D(vtk_file,"reflectivity","xray")

    # Add the clusters
    visit2D.add_clusters_2D(basename,"_cluster_",col_tables)

    # Add modes as labels
    modes_file=os.path.splitext(basename)[0]+"-clusters_modes.vtk"
    visit2D.add_modes(modes_file)

    # Get it all processed and stowed away
    DrawPlots()
    visit2D.save_window("tracked_",1)

    #
    # clean up
    #

    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(vtk_file)
    CloseDatabase(modes_file)
    visit2D.close_clusters_2D(basename,"_cluster_")
    return_code=call("rm *cluster_*.vtk", shell=True)

    # don't forget to increment run counter
    run_count = run_count + 1

#
# Create movies
#

print "Creating movies from slides"
return_code=call("convert -delay 50 source_*.png source.mpeg")
return_code=call("convert -delay 50 tracked_*.png tracked.mpeg")
return_code=call("convert -delay 50 untracked_*.png untracked.mpeg")

print "Done. Closing Visit."
Close();
