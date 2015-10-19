#!/usr/bin/python

MEANIE3D_HOME     = "P_M3D_HOME"
NETCDF_DIR        = "P_SOURCE_DIR"
PARAM_T           = "P_SCALE"
RESUME            = "P_RESUME"

try:
    import sys
    # Appending the module path is crucial
    sys.path.append(MEANIE3D_HOME+"/visualisation/modules")
    import os
    import meanie3D
    from subprocess import call
except ImportError as e:
    print 'Exception error is: %s' % e
    sys.exit()

#VAR_NAME="xband_oase_zh"
VAR_NAME="cband_oase_zh"

DETECT_PARAMS      = " -s "+PARAM_T
DETECT_PARAMS     += " --lower-thresholds "+VAR_NAME+"=20"
DETECT_PARAMS     += " --upper-thresholds "+VAR_NAME+"=65"
DETECT_PARAMS     += " -m 64"

CLUSTERING_PARAMS =  "-d z,y,x --vtk-dimensions x,y,z"
CLUSTERING_PARAMS += " --verbosity 1"
CLUSTERING_PARAMS += " --weight-function default"
CLUSTERING_PARAMS += " -v "+VAR_NAME
CLUSTERING_PARAMS += " " + DETECT_PARAMS

TRACKING_PARAMS =  " --verbosity 1 "
TRACKING_PARAMS += " -t "+VAR_NAME
TRACKING_PARAMS += " --wr=1.0 --ws=0.0 --wt=0.0"

meanie3D.run_tracking(NETCDF_DIR,".",CLUSTERING_PARAMS,TRACKING_PARAMS,bool(RESUME=="YES"),-1)
