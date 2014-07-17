#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "PARAM_SOURCE_DIR"
PARAM_T           = "PARAM_SCALE"
DYLD_PATH         = "PARAM_DL_PATH"

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

CLUSTERING_PARAMS =  " -s "+PARAM_T
CLUSTERING_PARAMS += " -v msevi_l15_ir_108 --upper-thresholds msevi_l15_ir_108=71.9493"
CLUSTERING_PARAMS += " -d y,x --vtk-dimensions x,y"
CLUSTERING_PARAMS += " --weight-function oase --wwf-lower-threshold 5"
CLUSTERING_PARAMS += " -m 10"
CLUSTERING_PARAMS += " --verbosity 1"

TRACKING_PARAMS =  " --vtk-dimensions x,y"
TRACKING_PARAMS += " -t msevi_l15_ir_108"
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"
TRACKING_PARAMS += " --verbosity 1 "

BIN_PREFIX = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_PATH+";"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

meanie3D.run_tracking(NETCDF_DIR,".",CLUSTERING_PARAMS,TRACKING_PARAMS,last_completed_run_count,-1)