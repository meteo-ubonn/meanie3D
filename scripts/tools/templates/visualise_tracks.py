#!/usr/bin/python

__author__ = "Juergen Simon"
__email__ = "juergen.simon@uni-bonn.de"
__version__ = "1.5.4"

# This is a template for a python script to be created with
# template variables replaced and then executed inside of
# Visit. The purpose is to visualise tracks from .vtk files
# produced by meanie3D-trackstats --write-center-tracks-as-vtk

# Parameters
_TRACKS_DIR =  "P_TRACKS_DIR"
_M3D_HOME    = "P_M3D_HOME"
_CONFIG_FILE = "P_CONFIGURATION_FILE"

# Import modules
import sys
sys.path.append(_M3D_HOME+"/scripts/python-modules")
import meanie3D
import meanie3D_visit

# Parse configuration data
configuration = meanie3D.load_configuration(_CONFIG_FILE)

# Add parsed parameters
configuration["meanie3d_home"] = _M3D_HOME
configuration["tracks_dir"] = _TRACKS_DIR

# run it
meanie3D_visit.plotTracks(configuration)