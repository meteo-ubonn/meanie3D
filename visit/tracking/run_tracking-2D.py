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
    import meanie3D
    from subprocess import call
except ImportError as e:
    print 'Exception error is: %s' % e
    sys.exit()

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# RADOLAN

VAR_NAME="RX"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds RX=46 -m 15"

CLUSTERING_PARAMS =  "-d y,x --vtk-dimensions x,y"
CLUSTERING_PARAMS += " --verbosity 1"
CLUSTERING_PARAMS += " --write-variables-as-vtk="+VAR_NAME
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " --write-clusters-as-vtk"
CLUSTERING_PARAMS += " --write-cluster-centers"

CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS

TRACKING_PARAMS = "--verbosity 1 --write-vtk"
TRACKING_PARAMS += " --vtk-dimensions x,y"
TRACKING_PARAMS += " -t "+VAR_NAME
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

meanie3D.run_tracking(NETCDF_DIR,".",CLUSTERING_PARAMS,TRACKING_PARAMS,last_completed_run_count,-1)
