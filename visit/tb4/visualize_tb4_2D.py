#!/usr/bin/python

# Parameters
NETCDF_FILE = "P_NETCDF_FILE"
CLUSTER_FILE = "P_CLUSTER_FILE"
VAR_NAME = "P_VAR_NAME"
M3D_HOME = "P_M3D_HOME"

# Import modules
import sys
sys.path.append(M3D_HOME+"/visit/modules")
import glob
import os
import time
import visit2D
import visitUtils
from subprocess import call

# Presumptuous

# TB4 files are usually satellite brightness temperatures.
# Cap at -20C for all intents

VAR_MAX=253.15

# Control

WITH_BACKGROUND_GRADIENT=False
WITH_TOPOGRAPHY=False
WITH_RIVERS_AND_BOUNDARIES=False
WITH_SOURCE_BACKROUND=False
WITH_DATETIME=False

CREATE_SOURCE_MOVIE=False
CREATE_CLUSTERS_MOVIE=False

# Conversion program params

CONVERSION_PARAMS  = "-t cluster "
CONVERSION_PARAMS += " -v "+VAR_NAME
CONVERSION_PARAMS += " --write-as-xml=false"
CONVERSION_PARAMS += " --extract-skin=false"
CONVERSION_PARAMS += " --vtk-dimensions lon,lat"

# binaries
DYLD_LIBRARY_PATH="/usr/local/lib:/usr/lib"
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"
print "Conversion Command:" + conversion_bin + " " + CONVERSION_PARAMS

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes

a = GetAnnotationAttributes()
a.axes2D.visible=1
a.axesArray.autoSetScaling=0
a.axes2D.xAxis.title.visible=1
a.axes2D.yAxis.title.visible=1
a.legendInfoFlag=1
a.databaseInfoFlag=0
a.userInfoFlag=0
a.timeInfoFlag=0
SetAnnotationAttributes(a)

# Add gray/black background gradient
if WITH_BACKGROUND_GRADIENT:
    visitUtils.add_background_gradient();

print "Cleaning up *.vt* *.png"
#return_code=call("rm -f *.vt* *.png", shell=True)

print "-- Creating colortables ---"
num_colors = visitUtils.create_cluster_colortable("cluster_colors")

# check if the files both exist
print "Visualzing file "+NETCDF_FILE+" and cluster file "+CLUSTER_FILE

if not os.path.exists(NETCDF_FILE):
    print "File "+NETCDF_FILE+" does not exist."
    exit(0)

if not os.path.exists(CLUSTER_FILE):
    print "File "+CLUSTER_FILE+" does not exist."
    exit(0)

# construct the cluster filename and find it
# in the cluster directory

netcdf_path,filename    = os.path.split(NETCDF_FILE);
basename                = os.path.splitext(filename)[0]

cluster_path,cluster_fn = os.path.split(CLUSTER_FILE)
cluster_basename        = os.path.split(cluster_fn)[0]
label_file              = cluster_basename+"-clusters_centers.vtk"

v = GetView2D()
v.windowCoords = (2,18.95,45,55.95)
v.viewportCoords = (0.2,0.95,0.15,0.95)
SetView2D(v)

OpenDatabase(NETCDF_FILE);

print "-- Plotting source data --"
start_time = time.time()

visit2D.add_pseudocolor(NETCDF_FILE,VAR_NAME,"hot_desaturated",1.0,1)
p = PseudocolorAttributes()
p.minFlag,p.maxFlag=1,1
p.max=VAR_MAX
p.invertColorTable=1
SetPlotOptions(p)

AddOperator("Threshold")
t = ThresholdAttributes();
t.upperBounds=(VAR_MAX)
SetOperatorOptions(t)

DrawPlots();
visitUtils.save_window(basename+"_"+VAR_NAME+"_",1)

DeleteAllPlots()
ClearWindow()

print "    done. (%.2f seconds)" % (time.time()-start_time)

print "-- Converting clusters to .vtr --"
start_time = time.time()
command=conversion_bin+" -f "+CLUSTER_FILE+" "+CONVERSION_PARAMS
print command
return_code = call( command, shell=True)
print "    done. (%.2f seconds)" % (time.time()-start_time)

print "-- Rendering cluster scene --"
start_time = time.time()

visit2D.add_pseudocolor(NETCDF_FILE,VAR_NAME,"hot_desaturated",0.1,1)
p = PseudocolorAttributes()
p.minFlag,p.maxFlag=1,1
p.max=VAR_MAX
p.invertColorTable=1
SetPlotOptions(p)

# Add the clusters
visit2D.add_clusters_with_colortable(cluster_basename,"_cluster_","cluster_colors",num_colors)

# Add modes as labels
visitUtils.add_labels(label_file,"geometrical_center")

# save as image
DrawPlots()

visitUtils.save_window(basename+"_"+VAR_NAME+"-clusters_",1)

print "    done. (%.2f seconds)" % (time.time()-start_time)

# clean up
DeleteAllPlots();
ClearWindow()
CloseDatabase(NETCDF_FILE)
CloseDatabase(label_file)
visit2D.close_topography()
visitUtils.close_pattern(basename+"*.vtr")
visitUtils.close_pattern(basename+"*.vtk")

#return_code=call("rm -f *.vt*", shell=True)
