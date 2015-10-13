#!/usr/bin/python

__author__ = "Juergen Simon"
__email__ = "juergen.simon@uni-bonn.de"
__version__ = "1.5.4"

# This is a template for a python script to be created with
# template variables replaced and then executed inside of
# Visit. The purpose is to visualise tracks from .vtk files
# produced by meanie3D-trackstats --write-center-tracks-as-vtk

# Parameters
_NETCDF_DIR     = "P_NETCDF_DIR"
_CLUSTER_DIR    = "P_CLUSTER_DIR"
_M3D_HOME       = "P_M3D_HOME"
_RESUME         = "P_RESUME"
_CONFIG_FILE    = "P_CONFIGURATION_FILE"

# Import modules
import sys
sys.path.append(_M3D_HOME)
import meanie3D.app.utils
import meanie3D.visualisation.clusters

# Parse configuration data
configuration = meanie3D.app.utils.load_configuration(_CONFIG_FILE)

# Add parsed parameters
configuration["meanie3d_home"] = _M3D_HOME
configuration["cluster_directory"] = _CLUSTER_DIR
configuration["source_directory"] = _NETCDF_DIR
configuration['resume'] = meanie3D.app.utils.strToBool(_RESUME)
configuration['config_file'] = _CONFIG_FILE

# run it
# pdb.set_trace()
meanie3D.visualisation.clusters.run(configuration)
quit()