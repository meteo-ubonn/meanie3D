#!/usr/bin/python

# Parameters
MEANIE3D_HOME     = "M3D_HOME"
NETCDF_DIR        = "SOURCE_DIR"
VAR_NAME          = "RADOLAN_VAR_NAME"

# Import modules
import sys
sys.path.append(MEANIE3D_HOME+"/visit/modules")
import glob
import os
import time
import visit2D
import visitUtils
from subprocess import call

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

print "SOURCE_DIR="+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME

RX_MIN=15
RX_MAX=75

# Set view and annotation attributes

print "Setting rendering attributes"
r = RenderingAttributes()
#r.geometryRepresentation = r.Points;
r.geometryRepresentation = r.Wireframe;
SetRenderingAttributes(r);

print ColorTableNames()

#print "Setting annotation attributes:"
#visit2D.set_annotations()

print "Cleaning up *.vtk *.vtr *.png *.gif"
return_code=call("rm -f *.vtk *.vtr *.png *.gif", shell=True)

# Setting 3D view parameters
print "Setting 2D view parameters"
visit2D.set_view_to_radolan();

# Add gray/black background gradient
#print "Setting background gradient"
#visitUtils.add_background_gradient();

#print "Creating colortables"
#visitUtils.create_topography_colortable()

# Set view and annotation attributes
a = GetAnnotationAttributes()
a.axes3D.visible=1
a.axes3D.autoSetScaling=0
a.userInfoFlag=0
a.timeInfoFlag=0
a.legendInfoFlag=1
a.databaseInfoFlag=1
a.axes3D.xAxis.title.visible=0
a.axes3D.yAxis.title.visible=0
a.axes3D.zAxis.title.visible=0
SetAnnotationAttributes(a)

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=sorted(glob.glob(netcdf_pattern))

run_count = 0

# Process the files one by one
for netcdf_file in netcdf_list:
    
    print "----------------------------------------------------------"
    print "Processing " + netcdf_file
    print "----------------------------------------------------------"
    
    #
    # Plot the source data in color
    #
    
    print "-- Adding map data --"
    visit2D.add_map_rivers("national")
    visit2D.add_map_borders("national")

    # plot data
    
    OpenDatabase(netcdf_file)
    AddPlot("Pseudocolor", VAR_NAME)
    p = PseudocolorAttributes()
    p.colorTableName = "hot_desaturated"
    p.legendFlag=1
    p.lightingFlag=1
    p.opacity=1.0
    
    if VAR_NAME == "RX":
        
        p.minFlag,p.maxFlag = 1,1
        p.min,p.max = RX_MIN,RX_MAX
        
        # threshold
        AddOperator("Threshold")
        t = ThresholdAttributes();
        t.lowerBounds=(RX_MIN)
        t.upperBounds=(RX_MAX)
        SetOperatorOptions(t)
    
    SetPlotOptions(p)

    # date/time
    visitUtils.add_datetime(netcdf_file)

    DrawPlots()
        
    visitUtils.save_window(VAR_NAME+"-",1)

    # clean up
    DeleteAllPlots();

    CloseDatabase(netcdf_file)

    # don't forget to increment run counter
    run_count = run_count + 1

    # memory leak fix
    if run_count % 100 == 0:
        CloseComputeEngine()

print "Done. Closing Visit."
exit()

# create loops
print "Creating animated gif ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"*.png radolan-"+VAR_NAME+".gif"
return_code=call(convert_cmd, shell=True)
print "Creating mpeg ..."
convert_cmd="/usr/local/bin/convert -limit memory 4GB -delay 50 -quality 100 "+VAR_NAME+"*.png radolan-"+VAR_NAME+".m4v"
return_code=call(convert_cmd, shell=True)

# clean up
print "Cleaning up ..."
return_code=call("mkdir images", shell=True)
return_code=call("mv *.png images", shell=True)
return_code=call("rm -f *.vt* *.py", shell=True)

print "Done."


