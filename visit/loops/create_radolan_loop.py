#!/usr/bin/python

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "SOURCE_DIR"
VAR_NAME          = "RADOLAN_VAR_NAME"  

# Appending the module path is crucial

sys.path.append(MEANIE3D_HOME+"/visit/modules")

import glob
import os
import sys
import time
import visit2D
import visitUtils
from subprocess import call

print [key for key in locals().keys()
       if isinstance(locals()[key], type(sys)) and not key.startswith('__')]

# TODO: find a more elegant way to resume
# if > 0 a previous run is resumed
last_completed_run_count = 0

# print parameters

print "Creating loop from files in directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes
a = GetAnnotationAttributes()
a.axes2D.visible=1
a.axes2D.autoSetScaling=0
a.userInfoFlag=0
a.timeInfoFlag=0
a.legendInfoFlag=1
a.databaseInfoFlag=1
SetAnnotationAttributes(a)

# Modify view parameters
v = GetView2D()
v.windowCoords = (-523.462,376.538,-4658.64,-3758.64)
v.viewportCoords = (0.2,0.95,0.15,0.95)
SetView2D(v)

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:
    
    print "----------------------------------------------------------"
    print "Processing " + netcdf_file
    print "----------------------------------------------------------"
    
    #
    # Plot the source data in color
    #
    
    visit2D.add_pseudocolor( netcdf_file, VAR_NAME, "hot_desaturated",1 )

    if VAR_NAME == "RX":
        cp=PseudocolorAttributes();
        cp.minFlag,cp.maxFlag = 1,1
        cp.min,cp.max = -32.5,55.0
        SetPlotOptions(cp)
        
    DrawPlots()
        
    visitUtils.save_window(VAR_NAME+"-",1)
        
    # clean up
    DeleteAllPlots();

    # don't forget to increment run counter
    run_count = run_count + 1

    # memory leak fix
    if run_count % 100 == 0:
        CloseComputeEngine()

print "Done. Closing Visit."
exit()

