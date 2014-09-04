#!/usr/bin/python

# Parameters
NETCDF_DIR = "P_NETCDF_DIR"
CLUSTER_DIR = "P_CLUSTER_DIR"
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

# General parameters
VAR_NAME = "RX"
VAR_MIN = 15;
VAR_MAX = 65;

# Control

WITH_BACKGROUND_GRADIENT=False
WITH_TOPOGRAPHY=False
WITH_RIVERS_AND_BOUNDARIES=True
WITH_SOURCE_BACKROUND=False
WITH_DATETIME=True

CREATE_SOURCE_MOVIE=False
CREATE_CLUSTERS_MOVIE=True

# Conversion program params

CONVERSION_PARAMS  = "-t cluster "
CONVERSION_PARAMS += " -v "+VAR_NAME
CONVERSION_PARAMS += " --write-as-xml=false"
CONVERSION_PARAMS += " --extract-skin=false"
CONVERSION_PARAMS += " --vtk-dimensions x,y"

# binaries

DYLD_LIBRARY_PATH="/usr/local/lib:/usr/lib"
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes

print "Setting annotation attributes:"
visit2D.set_annotations()

print "Cleaning up *.vt* tracking_*.png source_*.png"
return_code=call("rm -f *.vt* tracking_*.png source_*.png", shell=True)

# Set view to nationwide composite
visit2D.set_view_to_radolan();

print "-- Creating colortables ---"
num_colors = visitUtils.create_cluster_colortable("cluster_colors")

if WITH_TOPOGRAPHY:
    visitUtils.create_topography_colortable()

if WITH_BACKGROUND_GRADIENT:
    visitUtils.add_background_gradient();

# Glob the netcdf directory
netcdf_files = sorted(glob.glob(NETCDF_DIR+"/*.nc"));

print "Processing files in directory " + NETCDF_DIR
#print netcdf_files

# Keep track of number of images to allow
# forced re-set in time to circumvent the
# Visit memory leak
image_count=0

for netcdf_file in netcdf_files:

    # construct the cluster filename and find it
    # in the cluster directory
    
    netcdf_path,filename    = os.path.split(netcdf_file);
    basename                = os.path.splitext(filename)[0]

    cluster_file            = CLUSTER_DIR+"/"+basename+"-clusters.nc"
    label_file              = basename+"-clusters_centers.vtk"

    # check if the files both exist
    print "Visualzing file "+netcdf_file+" and cluster file "+cluster_file
    if not os.path.exists(cluster_file):
        print "Cluster file does not exist. Skipping."
        continue

    # plot source data?

    OpenDatabase(netcdf_file);

    if CREATE_SOURCE_MOVIE:

        if WITH_TOPOGRAPHY:
            print "-- Adding topography data --"
            visit2D.add_topography("national_topo_2D")

        if WITH_RIVERS_AND_BOUNDARIES:
            print "-- Adding map data --"
            visit2D.add_map_rivers("national")
            visit2D.add_map_borders("national")
            
        if WITH_DATETIME:
            print "-- Adding timestamp --"
            visitUtils.add_datetime(netcdf_file)

        print "-- Plotting source data --"
        start_time = time.time()

        visit2D.add_pseudocolor(netcdf_file,VAR_NAME,"hot_desaturated",1.0,1)
        p = PseudocolorAttributes()
        p.minFlag,p.maxFlag=1,1
        p.min,p.max=VAR_MIN,VAR_MAX
        SetPlotOptions(p)
        DrawPlots();

        AddOperator("Threshold")
        t = ThresholdAttributes();
        t.lowerBounds=(VAR_MIN)
        t.upperBounds=(VAR_MAX)
        SetOperatorOptions(t)

        DrawPlots()
        visitUtils.save_window("source_",1)
        image_count=image_count+1;

        DeleteAllPlots()
        ClearWindow()

        print "    done. (%.2f seconds)" % (time.time()-start_time)

    if CREATE_CLUSTERS_MOVIE:

        start_time = time.time()
        print "-- Converting clusters to .vtr --"

        # build the clustering command
        command=conversion_bin+" -f "+cluster_file+" "+CONVERSION_PARAMS
        print command
        return_code = call( command, shell=True)

        print "    done. (%.2f seconds)" % (time.time()-start_time)

        if WITH_TOPOGRAPHY:
            print "-- Adding topography data --"
            visit2D.add_topography("national_topo_2D")

        if WITH_RIVERS_AND_BOUNDARIES:
            print "-- Adding map data --"
            visit2D.add_map_rivers("national")
            visit2D.add_map_borders("national")

        if WITH_DATETIME:
            print "-- Adding timestamp --"
            visitUtils.add_datetime(netcdf_file)

        print "-- Rendering cluster scene --"
        start_time = time.time()
        
        if WITH_SOURCE_BACKROUND:

            # re-plot source data as canvas
            visit2D.add_pseudocolor(netcdf_file, VAR_NAME, "gray", 0.33, 0 )

            # threshold as before
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.lowerBounds=(VAR_MIN)
            t.upperBounds=(VAR_MAX)
            SetOperatorOptions(t)
        
        # Add the clusters
        basename = CLUSTER_DIR+"/"
        visit2D.add_clusters_with_colortable(basename,"_cluster_","cluster_colors",num_colors)

        # Add modes as labels
        visitUtils.add_labels(label_file,"geometrical_center")

        # save as image
        DrawPlots()

        visitUtils.save_window("tracking_",1)
        image_count=image_count+1;

        print "    done. (%.2f seconds)" % (time.time()-start_time)

    # clean up
    DeleteAllPlots();
    ClearWindow()
    CloseDatabase(netcdf_file)
    CloseDatabase(label_file)
    visit2D.close_topography()
    visitUtils.close_pattern(basename+"*.vtr")
    visitUtils.close_pattern(basename+"*.vtk")
    return_code=call("rm -f *.vt*", shell=True)

    # periodically kill computing engine to
    # work around the memory leak fix
    if image_count % 100 == 0:
        CloseComputeEngine()

# create loops
if CREATE_SOURCE_MOVIE:
    visitUtils.create_movie("source_","source.gif")
    visitUtils.create_movie("source_","source.m4v")

if  CREATE_CLUSTERS_MOVIE:
    visitUtils.create_movie("tracking_","tracking.gif")
    visitUtils.create_movie("tracking_","tracking.m4v")

# clean up
print "Cleaning up ..."
return_code=call("mkdir images", shell=True)
return_code=call("mv *.png images", shell=True)
return_code=call("mkdir movies", shell=True)
return_code=call("mv *.gif *.m4v movies", shell=True)
return_code=call("rm -f *.vt* visitlog.py", shell=True)

quit()