#!/usr/bin/python

# Parameters
SOURCE_DIR = "SOURCE_DIR_P"
VAR_NAME = "VARIABLE_P"
COLOR_TABLE = "COLOR_TABLE_P"
INVERT_COLOR_TABLE = "COLOR_TABLE_INVERT_P"
M3D_HOME = "MEANIE3D_HOME_P"

# Import modules
import sys
sys.path.append(M3D_HOME+"/visit/modules")
import glob
import os
import time
import visit2D
import visitUtils
from subprocess import call

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

print "SOURCE_DIR="+SOURCE_DIR
print "VAR_NAME="+VAR_NAME
print "COLOR_TABLE="+COLOR_TABLE
print "INVERT_COLOR_TABLE="+INVERT_COLOR_TABLE
print "M3D_HOME="+M3D_HOME

# Set view and annotation attributes

print "Setting annotation attributes:"
visit2D.set_annotations()

#print "Cleaning up *.vtk *.vtr"
#return_code=call("rm -f *.vtk *.vtr", shell=True)

#print "Cleaning up "+VAR_NAME+"*.png"
#return_code=call("rm -f "+VAR_NAME+"*.png", shell=True)

# Setting 3D view parameters
print "Setting 2D view parameters"
visit2D.set_view_to_radolan();

# Add gray/black background gradient
print "Setting background gradient"
visitUtils.add_background_gradient();

print "Creating colortables"
visitUtils.create_topography_colortable()

# Glob the netcdf directory
print "Processing files in directory " + SOURCE_DIR
netcdf_files = glob.glob(SOURCE_DIR+"/*.nc");

for netcdf_file in netcdf_files:
    
    print "Proessing " + netcdf_file
    
    # add 2D topograpy
    visit2D.add_mapstuff("national")
    
    # now plot the data
    OpenDatabase(netcdf_file)
    
    # Plot source
    visit2D.add_pseudocolor(netcdf_file,VAR_NAME,COLOR_TABLE,1,1)

    if INVERT_COLOR_TABLE == 1:
        p = PseudocolorAttributes()
        p.invertColorTable = 1
        SetPlotOptions(p)

    # date/time
    visitUtils.add_datetime(netcdf_file)

    DrawPlots()
    CloseDatabase(netcdf_file)
    visit2D.close_mapstuff();
    visitUtils.save_window(VAR_NAME+"_",1)
    DeleteAllPlots()
    ClearWindow()

print "Creating movie ..."
return_code=call("mkdir "+VAR_NAME, shell=True)
return_code=call("mkdir "+VAR_NAME+"/images", shell=True)

convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"_*.png "+VAR_NAME+".gif"
return_code=call(convert_cmd, shell=True)

convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"_*.png "+VAR_NAME+".m4v"
return_code=call(convert_cmd, shell=True)

return_code=call("mv "+VAR_NAME+"*.png "+VAR_NAME+"/images", shell=True)
return_code=call("mv "+VAR_NAME+"*.gif "+VAR_NAME, shell=True)
return_code=call("mv "+VAR_NAME+"*.m4v "+VAR_NAME, shell=True)

print "Done."
exit(0)