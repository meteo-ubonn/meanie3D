#!/usr/bin/python

# Parameters
SOURCE_DIR = "SOURCE_DIR_P"
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
print "M3D_HOME="+M3D_HOME

# Set view and annotation attributes

print "Setting rendering attributes"
r = RenderingAttributes()
print r
r.geometryRepresentation = r.Points;
print r
SetRenderingAttributes(r);

print ColorTableNames()

print "Setting annotation attributes:"
visit2D.set_annotations()

print "Cleaning up *.vtk *.vtr *.png *.gif"
return_code=call("rm -f *.vtk *.vtr *.png *.gif", shell=True)

# Setting 3D view parameters
print "Setting 2D view parameters"
visit2D.set_view_to_radolan();

# Add gray/black background gradient
print "Setting background gradient"
visitUtils.add_background_gradient();

#print "Creating colortables"
#visitUtils.create_topography_colortable()

# Glob the netcdf directory
print "Processing files in directory " + SOURCE_DIR
netcdf_files = glob.glob(SOURCE_DIR+"/*.nc");

for netcdf_file in netcdf_files:
    
    print "Proessing " + netcdf_file
    
    # add 2D topograpy
    visit2D.add_mapstuff("national")
    
    # now plot the data
    OpenDatabase(netcdf_file)
    
    #
    # Satellite
    #
    
    # Plot source
    visit2D.add_pseudocolor(netcdf_file,"msevi_l15_ir_108","bluehot",1,1)
    
    # threshold
    AddOperator("Threshold")
    t = ThresholdAttributes();
    t.upperBounds=(48.57)
    SetOperatorOptions(t)

    # Invert and logarithmic scaling
    p = PseudocolorAttributes()
    p.colorTableName="bluehot"
    p.invertColorTable=1
    p.scaling = p.Skew;
    SetPlotOptions(p)
    
    #
    # Radar
    #
    
    # Plot source
    visit2D.add_pseudocolor(netcdf_file,"cband_radolan_rx","hot_desaturated",1,1)
    
    # threshold
    AddOperator("Threshold")
    t = ThresholdAttributes();
    t.lowerBounds=(70)
    t.upperBounds=(150)
    SetOperatorOptions(t)
    
    #
    # Lightning
    #
    
    # Plot source
    visit2D.add_pseudocolor(netcdf_file,"linet_oase_tl","difference",1,1)
    
    # threshold
    AddOperator("Threshold")
    t = ThresholdAttributes();
    t.lowerBounds=(1)
    SetOperatorOptions(t)

    # Invert and logarithmic scaling
    p = PseudocolorAttributes()
    p.colorTableName = "difference"
    p.pointSizePixels = 5
    p.invertColorTable = 1
    SetPlotOptions(p)
    

    # date/time
    visitUtils.add_datetime(netcdf_file)

    DrawPlots()
    visitUtils.save_window("oase_composite_",1)
    
    CloseDatabase(netcdf_file)
    visit2D.close_mapstuff();
    DeleteAllPlots()
    ClearWindow()


# create loops
print "Creating animated gif ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 oase_composite_*.png oase_composite.gif"
return_code=call(convert_cmd, shell=True)
print "Creating mpeg ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 oase_composite_*.png oase_composite.m4v"
return_code=call(convert_cmd, shell=True)

# clean up
print "Cleaning up ..."
return_code=call("mkdir composite", shell=True)
return_code=call("mkdir composite/images", shell=True)
return_code=call("mv *.png composite/images", shell=True)
return_code=call("mv *.gif composite", shell=True)
return_code=call("mv *.m4v composite", shell=True)

print "Done."
