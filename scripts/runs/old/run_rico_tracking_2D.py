#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "PARAM_SOURCE_DIR"
PARAM_T           = "PARAM_SCALE"

try:
    import sys
    # Appending the module path is crucial
    sys.path.append(MEANIE3D_HOME+"/visit/modules")
    import os
    from meanie3D import meanie3D
    from subprocess import call
except ImportError as e:
    print 'Exception error is: %s' % e
    sys.exit()

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# RADOLAN

VAR_NAME="lwp"

#DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS      = " -r "+PARAM_T+","+PARAM_T+",10"

CLUSTERING_PARAMS =  "-d yt,xt --vtk-dimensions xt,yt"
CLUSTERING_PARAMS += " --verbosity 2"
#CLUSTERING_PARAMS += " --write-variables-as-vtk="+VAR_NAME
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " --lower-thresholds lwp=0.1"
CLUSTERING_PARAMS += " --wwf-lower-threshold 0.0"
CLUSTERING_PARAMS += " --min-cluster-size 10"
#CLUSTERING_PARAMS += " --coalesce-with-strongest-neighbour"

CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS

#TRACKING_PARAMS = "--verbosity 1 --write-vtk"
TRACKING_PARAMS = "--verbosity 3"
TRACKING_PARAMS += " --vtk-dimensions xt,yt"
TRACKING_PARAMS += " -t "+VAR_NAME
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

for time_index in xrange(0, 100):

    CP = CLUSTERING_PARAMS + " --time-index " + str(time_index)

    meanie3D.run_tracking(NETCDF_DIR,
                          ".",
                          CP,
                          TRACKING_PARAMS,
                          last_completed_run_count,
                          time_index)
