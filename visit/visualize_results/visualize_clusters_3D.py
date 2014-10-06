#!/usr/bin/python

# Parameters
NETCDF_DIR = "P_NETCDF_DIR"
CLUSTER_DIR = "P_CLUSTER_DIR"
M3D_HOME = "P_M3D_HOME"
RESUME = "P_RESUME"

# Import modules
import sys
sys.path.append(M3D_HOME+"/visit/modules")
import glob
import os
import time
import visit3D
import visitUtils
from subprocess import call

<<<<<<< HEAD
# Silent
SuppressMessages(True)
SuppressQueryOutputOn()
=======
# General parameters
VAR_NAME="cband_oase_zh"
VAR_MIN=20
VAR_MAX=65
>>>>>>> be73d22ddb81adf3bdb267305032c2f07616c71e

# Control

WITH_BACKGROUND_GRADIENT=False
WITH_TOPOGRAPHY=False
WITH_RIVERS_AND_BOUNDARIES=True
WITH_SOURCE_BACKROUND=False
WITH_DATETIME=True
CREATE_SOURCE_MOVIE=True
CREATE_CLUSTERS_MOVIE=True

# Stretching data and objects for better 3D impression
SCALE_FACTOR_Z=3.0

# 'local'    covers cologne/bonn/juelich area.
# 'national' covers Germany area
GRID_EXTENT="national"

# Conversion program params

# variables
variables=["cband_oase_zh"]

# thresholds
lower_thresholds=[35.0]
upper_thresholds=[65.0]

# legend scale
VAR_MIN=[0.0]
VAR_MAX=[65.0]

colortables=["hot_desaturated"]
colortables_invert_flags=[0]
opacity=[0.33]

CONVERSION_PARAMS  = "-t cluster "
CONVERSION_PARAMS += " -v "+variables[0]
CONVERSION_PARAMS += " --write-as-xml=false"
CONVERSION_PARAMS += " --extract-skin=false"
CONVERSION_PARAMS += " --vtk-dimensions x,y,z"

configuration = {
    'NETCDF_DIR' : NETCDF_DIR,
    'CLUSTER_DIR' : CLUSTER_DIR,
    'M3D_HOME' : M3D_HOME,
    'RESUME' : bool(RESUME=="YES"),
    'WITH_BACKGROUND_GRADIENT' : WITH_BACKGROUND_GRADIENT,
    'WITH_TOPOGRAPHY' : WITH_TOPOGRAPHY,
    'WITH_RIVERS_AND_BOUNDARIES' : WITH_RIVERS_AND_BOUNDARIES,
    'WITH_SOURCE_BACKROUND' : WITH_SOURCE_BACKROUND,
    'WITH_DATETIME' : WITH_DATETIME,
    'CREATE_SOURCE_MOVIE' : CREATE_SOURCE_MOVIE,
    'CREATE_CLUSTERS_MOVIE' : CREATE_CLUSTERS_MOVIE,
    'SCALE_FACTOR_Z' : SCALE_FACTOR_Z,
    'GRID_EXTENT' : GRID_EXTENT,
    'CONVERSION_PARAMS' : CONVERSION_PARAMS,
    'VARIABLES' : variables,
    'LOWER_TRESHOLDS' : lower_thresholds,
    'UPPER_TRESHOLDS' : upper_thresholds,
    'VAR_MIN' : VAR_MIN,
    'VAR_MAX' : VAR_MAX,
    'COLORTABLES' : colortables,
    'COLORTABLES_INVERT_FLAGS' : colortables_invert_flags,
    'OPACITY' : opacity
}

visit3D.create_multiperspective(configuration)

quit()