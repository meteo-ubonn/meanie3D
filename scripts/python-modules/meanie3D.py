#!/usr/bin/python

# -------------------------------------------------------------------
# Filename: meanie3D.py
#
# Contains some utility functions and the code that runs complete
# directories through the clustering/tracking process
# -------------------------------------------------------------------

import glob
import os
import os.path
import string
import shutil
import time
import json
from subprocess import call

# -------------------------------------------------------------------
# Define some executables to be called from python
# -------------------------------------------------------------------

BIN_PREFIX = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/usr/local/lib;"
detection_bin = BIN_PREFIX + "meanie3D-detect"
tracking_bin  = BIN_PREFIX + "meanie3D-track"
trackplot_bin = BIN_PREFIX + "meanie3D-trackplot"

# -------------------------------------------------------------------
# @return current module version
# -------------------------------------------------------------------

def get_version():
    return "1.5.0"

# -------------------------------------------------------------------
# parses a JSON configuration file
# @param filename
# @return configuration dictionary
# -------------------------------------------------------------------
def load_configuration(filename):
    json_data=open(filename)
    data = json.load(json_data)
    json_data.close()
    return data;

# -------------------------------------------------------------------
# Counts the number of netcdf files in the given
# directory
# @param directory
# @return number of netcdf files
# -------------------------------------------------------------------
def number_of_netcdf_files(source_dir):
    netcdf_pattern = source_dir + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    return len(netcdf_list)

# -------------------------------------------------------------------
# Creates an output filename based on given filename by
# appending -<slicenum>.nc at the end.
# @param basic filename
# @param slice num 
# @return filename-1.nc
# -------------------------------------------------------------------
def numbered_filename(filename,index):
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0]+"-"+str(index)+".nc"
    
# -------------------------------------------------------------------
# Deletes the directories 'log' and 'netcdf' underneath
# base path. Removes previous ones if they do exist
#
# @param base path
# -------------------------------------------------------------------
def create_ouput_directories(base_path):

    # base path

    if not os.path.exists(base_path):
        os.makedirs(base_path)

    # logs

    log_dir = base_path+"/log"
    if os.path.exists(log_dir):
        shutil.rmtree(log_dir)
    os.makedirs(log_dir)
    
    # netcdf results

    netcdf_dir = base_path+"/netcdf"
    if os.path.exists(netcdf_dir):
        shutil.rmtree(netcdf_dir)
    os.makedirs(netcdf_dir)

    return

# -------------------------------------------------------------------
# Runs a batch of files through the clustering and tracking.
#
# @param configuration with the following keys:
# 
#    DESCRIPTION : a description of the configuration
#    NETCDF_DIR : directory containing files to process
#    OUTPUT_DIR : directory to write results to
#    M3D_HOME   : ${MEANIE3D_HOME} environment variable
#    RESUME     : if True then pick off where you left off, if False
#                 all previous results are deleted and processing
#                 starts from scratch
#    CLUSTERING_PARAMS : parameters for meanie3D-detect
#    TRACKING_PARAMS: parameters for meanie3D-track
#    SCALE : scale parameter (optional)
#    USE_PREVIOUS: True: using previous results to enhance tracking
#                  Absent or False: no use of previous results.
#    USE_CI_SCORE : Use the CI-Score algorithm
#
# @param time index the index in time dimension to read and process.
#
# -------------------------------------------------------------------
def run_tracking(config,time_index):

    print "---------------------------------------------------"
    print "Configuration: " + config['DESCRIPTION']
    print "---------------------------------------------------"

    resume_at_index = 0;

    # In case a scale T parameter is given, the output dir
    # is scale<T>. Otherwise it's 'clustering'

    output_dir = config['OUTPUT_DIR']
    if (config.get('SCALE') == "None"):
        output_dir = output_dir + "/clustering"
    else:
        output_dir = output_dir + "/scale"+str(config['SCALE'])

    print "Writing output to " + output_dir

    # Resume?

    if config['RESUME'] == False:

        # consider time index. Even if not resuming, the 
        # output directories should only be created at
        # the first time step

        if time_index <= 0: 
            print "Removing results from previous runs"
            create_ouput_directories(output_dir)

    else:

        resume_at_index = number_of_netcdf_files(output_dir+"/netcdf")
        print "Resuming at index " + str(resume_at_index)

    # Get a list of the files we need to process

    use_ci_score = False
    if "USE_CI_SCORE" in config.keys():
        use_ci_score = config["USE_CI_SCORE"]

    netcdf_pattern = config['NETCDF_DIR'] + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    last_cluster_file=""
    run_count = 0

    # Process the files one by one
    
    for netcdf_file in netcdf_list:
        
        basename = os.path.basename(netcdf_file)
        
        cluster_file= ""
        
        if time_index < 0:
            cluster_file= output_dir+"/netcdf/" + os.path.splitext(basename)[0] + "-clusters.nc"
        else:
            cluster_file= output_dir+"/netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index) + ".nc"
            last_cluster_file = output_dir+"/netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index-1) + ".nc"
                
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

        # build the clustering command

        command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + config['CLUSTERING_PARAMS'] 

        # amend for time index and make proper logfile name

        if time_index < 0:
            logfile = output_dir+"/log/clustering_" + str(run_count)+".log"
        else:
            command = command + " -t " + str(time_index)
            logfile = output_dir+"/log/clustering_" + str(time_index)+".log"
            
        # scale?

        if not config.get('SCALE') == "None":
            scale_param = " -s " + config['SCALE']
            command += scale_param
        
        # use previous result to enhance current?

        if ((run_count > 0) or (time_index > 0)) and config['USE_PREVIOUS']:
            command += " -p " + last_cluster_file

        # add ci-comparison-file if applicable

        if run_count >= 3 and use_ci_score:
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
        
        #
        # Tracking
        #
        
        # if we have a previous scan, run the tracking command
        
        if (run_count > 0) or (time_index > 0):
            
            if time_index < 0:
                logfile = output_dir+"/log/tracking_" + str(run_count)+".log"
            else:
                logfile = output_dir+"/log/tracking_" + str(time_index)+".log"

            print "-- Tracking --"
            command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + config['TRACKING_PARAMS']
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
