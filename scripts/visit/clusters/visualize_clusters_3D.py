#!/usr/bin/python

# Parameters
NETCDF_DIR = "P_NETCDF_DIR"
CLUSTER_DIR = "P_CLUSTER_DIR"
M3D_HOME = "P_M3D_HOME"
RESUME = "P_RESUME"
CONFIG_FILE = "P_CONFIGURATION_FILE"

# Import modules
import sys
sys.path.append(M3D_HOME+"/scripts/python-modules")
from meanie3D import visit3D
import meanie3D

# Parse configuration data
configuration = meanie3D.load_configuration(CONFIG_FILE);

# Add parsed parameters
configuration["NETCDF_DIR"] = NETCDF_DIR
configuration["CLUSTER_DIR"] = CLUSTER_DIR
configuration["M3D_HOME"] = M3D_HOME
configuration["RESUME"] = bool(RESUME=="YES")

# run it
visit3D.visualization(configuration)

quit()