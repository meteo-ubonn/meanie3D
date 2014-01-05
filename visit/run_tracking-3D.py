#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "SOURCE_DIR"
PARAM_T           = "SCALE"

# Appending the module path is crucial

try:
    import sys
    sys.path.append(MEANIE3D_HOME+"/visit/modules")
    
    import glob
    import os
    import time
    import shutil
    
    from subprocess import call

except ImportError as e:
    print 'Exception error is: %s' % e
    sys.exit()

# Set up directory structure for results

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# if not resuming, create direcories
if last_completed_run_count == 0:
    
    # logs
    if os.path.exists('log'):
        shutil.rmtree('log')
    os.makedirs('log')

    # results
    if os.path.exists('netcdf'):
        shutil.rmtree('netcdf')
    os.makedirs('netcdf')

# RADOLAN

VAR_NAME="zh"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds "+VAR_NAME+"=35"
DETECT_PARAMS     += " -m 10"

CLUSTERING_PARAMS =  "-d z,y,x --vtk-dimensions x,y,z"
CLUSTERING_PARAMS += " --verbosity 1"
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS

TRACKING_PARAMS =  " --verbosity 1 "
TRACKING_PARAMS += " -t "+VAR_NAME
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

toggled_maintain=False

# delete previous results
if last_completed_run_count == 0:
    print "Cleaning up *.vtk *.nc *.png *.log"
    return_code=call("rm -f *.nc *.vtk *.log *.png", shell=True)

# binaries

detection_bin = "meanie3D-detect"
tracking_bin  = "meanie3D-track"
trackplot_bin = "meanie3D-trackplot"

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
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
    
    print "----------------------------------------------------------"
    print "Processing " + netcdf_file
    print "----------------------------------------------------------"
    
    print "-- Clustering --"

    #
    # Cluster
    #
    
    # build the clustering command
    command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + CLUSTERING_PARAMS

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
        command =tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + TRACKING_PARAMS
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
exit()
