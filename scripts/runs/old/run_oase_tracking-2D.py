#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "PARAM_SOURCE_DIR"
PARAM_T           = "PARAM_SCALE"
DYLD_PATH         = "PARAM_DL_PATH"
RESUME            = "P_RESUME"

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

CLUSTERING_PARAMS =  " -s "+PARAM_T
#CLUSTERING_PARAMS += " -v msevi_l15_ir_108 --upper-thresholds msevi_l15_ir_108=71.9493"
#CLUSTERING_PARAMS += " --weight-function oase --wwf-lower-threshold 5"
CLUSTERING_PARAMS += " -v msevi_l15_ir_108 --upper-thresholds msevi_l15_ir_108=28.27"
CLUSTERING_PARAMS += " --weight-function inverse --wwf-lower-threshold 0"

CLUSTERING_PARAMS += " -d y,x --vtk-dimensions x,y"
CLUSTERING_PARAMS += " -m 10"
CLUSTERING_PARAMS += " --verbosity 1"

TRACKING_PARAMS =  " --vtk-dimensions x,y"
TRACKING_PARAMS += " -t msevi_l15_ir_108"
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"
TRACKING_PARAMS += " --verbosity 1 "

meanie3D.run_tracking(NETCDF_DIR,".",CLUSTERING_PARAMS,TRACKING_PARAMS,bool(RESUME=="YES"),-1)