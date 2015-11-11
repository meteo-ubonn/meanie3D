#!/usr/bin/python

__author__ = "Juergen Simon"
__email__ = "juergen.simon@uni-bonn.de"
__version__ = "1.5.4"

# This is a template for a python script to be created with
# template variables replaced and then executed inside of
# Visit. The purpose is to visualise tracks from .vtk files
# produced by meanie3D-trackstats --write-center-tracks-as-vtk

# Parameters
_TRACKS_DIR     = "P_TRACKS_DIR"
_M3D_HOME       = "P_M3D_HOME"
_CONFIG_FILE    = "P_CONFIGURATION_FILE"
_RESUME         = "P_RESUME"

# Import modules
import os
import pdb
import sys

sys.path.append(_M3D_HOME)
import meanie3D.app.utils
import meanie3D.visualisation.tracks

# Parse configuration data
configuration = meanie3D.app.utils.load_configuration(_CONFIG_FILE)

# Add parsed parameters
configuration["meanie3d_home"] = _M3D_HOME
configuration["tracks_dir"] = _TRACKS_DIR
configuration['resume'] = meanie3D.app.utils.strToBool(_RESUME)

# run it
# pdb.set_trace()
meanie3D.visualisation.tracks.run(configuration)
quit()