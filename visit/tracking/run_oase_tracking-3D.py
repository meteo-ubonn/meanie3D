#!/usr/bin/python

# python script template for tracking

# TODO: parse command line params to this script instead of templating
# TODO: Extend this script to handle multiple scales at once

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "SOURCE_DIR"
PARAM_T           = "SCALE"

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

# RADOLAN

VAR_NAME="xband_oase_zh"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds "+VAR_NAME+"=35"
DETECT_PARAMS     += " --upper-thresholds "+VAR_NAME+"=75"
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

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

meanie3D.run_tracking(NETCDF_DIR,".",CLUSTERING_PARAMS,TRACKING_PARAMS,last_completed_run_count,-1)

