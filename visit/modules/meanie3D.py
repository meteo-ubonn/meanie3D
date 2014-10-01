#!/usr/bin/python
# Filename: meanie3D.py

version = 'v1.4'

# This module contains some core functionality that
# is re-used all over

import glob
import os
import os.path
import string
import shutil
import time
from subprocess import call

# TODO: *sigh*
BIN_PREFIX = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/usr/local/lib;"

detection_bin = BIN_PREFIX + "meanie3D-detect"
tracking_bin  = BIN_PREFIX + "meanie3D-track"
trackplot_bin = BIN_PREFIX + "meanie3D-trackplot"

# Counts the number of netcdf files in the given
# directory
def number_of_netcdf_files(source_directory):
    netcdf_pattern = source_directory + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    return len(netcdf_list)

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
def run_tracking(source_directory,
                 output_directory,
                 clustering_params,
                 tracking_params,
                 resume,
                 time_index):
    
    if resume:
        
        # find the index to resume from by counting the
        # files in the result directory
        
        resume_at_index = number_of_netcdf_files(output_directory+"/netcdf")
    
    else:

        # delete previous results (if not in a time series)
        
        resume_at_index = 0;

        if time_index < 0:
            
            # if not resuming, create direcories
            print "Creating output directory structure"
            create_per_scale_directories(output_directory)

            print "Cleaning up *.vtk *.nc *.png *.log"
            return_code=call("rm -f *.nc *.vtk *.log *.png", shell=True)
        
        elif time_index == 0:
            
            # at first time step, create directories
            print "Creating output directory structure"
            create_per_scale_directories(output_directory)


    # Get a list of the files we need to process
    netcdf_pattern = source_directory + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))

    last_cluster_file=""
    run_count = 0

    # Process the files one by one
    
    for netcdf_file in netcdf_list:
        
        basename = os.path.basename(netcdf_file)
        
        cluster_file= ""
        
        if time_index < 0:
            cluster_file= "./netcdf/" + os.path.splitext(basename)[0] + "-clusters.nc"
        else:
            cluster_file= "./netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index) + ".nc"
            last_cluster_file = "./netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index-1) + ".nc"
        
        
        # if there is a resume counter, keep skipping
        # until the count is right
        if (resume_at_index > 0) and (run_count < resume_at_index):
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

        if time_index < 0:
            logfile = "./log/clustering_" + str(run_count)+".log"
        else:
            logfile = "./log/clustering_" + str(time_index)+".log"
        
        # build the clustering command
        command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + clustering_params
        
        # use previous result to enhance current
        if (run_count > 0) or (time_index > 0):
            command = command + " -p " + last_cluster_file
        
        # complete command with directing output to logfile
        command = command + " > " + logfile
        
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
        
        if (run_count > 0) or (time_index > 0):
            
            if time_index < 0:
                logfile = "./log/tracking_" + str(run_count)+".log"
            else:
                logfile = "./log/tracking_" + str(time_index)+".log"

            print "-- Tracking --"
            command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + tracking_params
            command = command + " > " + logfile
            
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


# Runs a complete directory through the tracking
# using the meanie3d-detect-ci algorithm
#
# @param source directory
# @param output directory
# @param flags for detection program
# @param flags for tracking program
# @param count of last completed file (for resuming)
#
def run_ci_tracking(source_directory,
                    output_directory,
                    clustering_params,
                    tracking_params,
                    resume,
                    time_index):
    
    if resume:
        
        # find the index to resume from by counting the
        # files in the result directory
        
        resume_at_index = number_of_netcdf_files(output_directory+"/netcdf")
    
    else:
        
        # delete previous results (if not in a time series)
        
        resume_at_index = 0;
        
        if time_index < 0:
            
            # if not resuming, create direcories
            print "Creating output directory structure"
            create_per_scale_directories(output_directory)
            
            print "Cleaning up *.vtk *.nc *.png *.log"
            return_code=call("rm -f *.nc *.vtk *.log *.png", shell=True)
        
        elif time_index == 0:
            
            # at first time step, create directories
            print "Creating output directory structure"
            create_per_scale_directories(output_directory)
    
    
    # Get a list of the files we need to process
    netcdf_pattern = source_directory + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    
    last_cluster_file=""
    run_count = 0
    
    # Process the files one by one
    
    for netcdf_file in netcdf_list:
        
        basename = os.path.basename(netcdf_file)
        
        cluster_file= ""
        
        if time_index < 0:
            cluster_file= "./netcdf/" + os.path.splitext(basename)[0] + "-clusters.nc"
        else:
            cluster_file= "./netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index) + ".nc"
            last_cluster_file = "./netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index-1) + ".nc"
        
        
        # if there is a resume counter, keep skipping
        # until the count is right
        if (resume_at_index > 0) and (run_count < resume_at_index):
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
        
        if time_index < 0:
            logfile = "./log/clustering_" + str(run_count)+".log"
        else:
            logfile = "./log/clustering_" + str(time_index)+".log"
        
        # build the clustering command
        command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + clustering_params
        
        # use previous result to enhance current
        if (run_count > 0) or (time_index > 0):
            command = command + " -p " + last_cluster_file
        
        # add ci-comparison-file if applicable
        if run_count >= 3:
            command += " --ci-comparison-file " + netcdf_list[run_count-3]
            proto_file = os.path.splitext(os.path.basename(netcdf_list[run_count-3]))[0] + "-protoclusters.nc"
            command += " --ci-comparison-protocluster-file " + proto_file;
        
        # complete command with directing output to logfile
        command = command + " > " + logfile
        
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
        
        if (run_count > 0) or (time_index > 0):
            
            if time_index < 0:
                logfile = "./log/tracking_" + str(run_count)+".log"
            else:
                logfile = "./log/tracking_" + str(time_index)+".log"
            
            print "-- Tracking --"
            command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + tracking_params
            command = command + " > " + logfile
            
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