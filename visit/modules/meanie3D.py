#!/usr/bin/python
# Filename: meanie3D.py

version = 'v1.3'

# This module contains some core functionality that
# is re-used all over

import glob
import os
import os.path
import string
import shutil
import time
from subprocess import call

# binaries
detection_bin = "meanie3D-detect"
tracking_bin  = "meanie3D-track"
trackplot_bin = "meanie3D-trackplot"

# Deletes the directories 'log' and 'netcdf' underneath
# base path. Removes previous ones if they do exist
#
# @param base path
#
def create_per_scale_directories(base_path):

    # logs
    log_dir = base_path+"/log"
    if os.path.exists(log_dir):
        shutil.rmtree(log_dir)
    os.makedirs(log_dir)
    
    # results
    netcdf_dir = base_path+"/netcdf"
    if os.path.exists(netcdf_dir):
        shutil.rmtree(netcdf_dir)
    os.makedirs(netcdf_dir)

    return

# Runs a complete directory through the tracking.
#
# @param source directory
# @param output directory
# @param flags for detection program
# @param flags for tracking program
# @param count of last completed file (for resuming)
#
def run_tracking(source_directory,output_directory,clustering_params,tracking_params,last_completed_run_count):

    # TODO: find a more elegant way to resume
    # if > 0 a previous run is resumed
    last_completed_run_count = 0

    if last_completed_run_count == 0:
        
        # if not resuming, create direcories
        print "Creating output directory structure"
        create_per_scale_directories(output_directory)

        # delete previous results
        print "Cleaning up *.vtk *.nc *.png *.log"
        return_code=call("rm -f *.nc *.vtk *.log *.png", shell=True)

    # Get a list of the files we need to process
    netcdf_pattern = source_directory + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    last_cluster_file=""

    run_count = 0

    # Process the files one by one
    for netcdf_file in netcdf_list:
        
        basename = os.path.basename(netcdf_file)
        cluster_file= "./netcdf/" + os.path.splitext(basename)[0] + "-clusters.nc"
        
        # if there is a resume counter, keep skipping
        # until the count is right
        if (last_completed_run_count > 0) and (run_count <= last_completed_run_count):
            last_cluster_file = cluster_file
            run_count = run_count + 1
            continue
        
        print "-------------------------------------------------------------------------------------"
        print "Processing " + netcdf_file
        print "-------------------------------------------------------------------------------------"

        print "-- Clustering --"
        
        #
        # Cluster
        #
        
        # build the clustering command
        command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + clustering_params
        
        # use previous result to enhance current
        if run_count > 0:
            command = command + " -p " + last_cluster_file
        
        command = command + " > ./log/clustering_" + str(run_count)+".log"
        
        # execute
        print command
        start_time = time.time()
        return_code = call( command, shell=True)
        print "    done. (%.2f seconds)" % (time.time()-start_time)
        
        # clean vtk/vtr files
        return_code=call("rm -f *cluster*_*.vt*", shell=True)
        
        #
        # Tracking
        #
        
        # if we have a previous scan, run the tracking command
        
        if run_count > 0:
            
            print "-- Tracking --"
            command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + tracking_params
            command = command + " > ./log/tracking_" + str(run_count)+".log"
            
            # execute
            print command
            start_time = time.time()
            return_code = call( command, shell=True)
            print "    done. (%.2f seconds)" % (time.time()-start_time)
        
        # keep track
        last_cluster_file=cluster_file
        
        print "Cleaning up *.vt*"
        return_code=call("rm -f *.vt*", shell=True)
        
        # don't forget to increment run counter
        run_count = (run_count + 1)

    print "Done."
    return

