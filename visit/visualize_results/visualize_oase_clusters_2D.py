#!/usr/bin/python

# Parameters
NETCDF_DIR = "P_NETCDF_DIR"
CLUSTER_DIR = "P_CLUSTER_DIR"
M3D_HOME = "P_M3D_HOME"
OVERLAY_VALUE = "P_OVERLAY"
RESUME_AT_FILE = "P_RESUME_AT_FILE"

# Import modules
import sys
sys.path.append(M3D_HOME+"/visit/modules")
import glob
import os
import time
#import visit2D
import visit2D
import visitUtils
from subprocess import call

#print [key for key in locals().keys()
#       if isinstance(locals()[key], type(sys)) and not key.startswith('__')]

OVERLAY=bool(OVERLAY_VALUE == "YES")
CREATE_SOURCE_MOVIE=True
CREATE_CLUSTERS_MOVIE=True
OVERLAY_CLUSTERS=True

WITH_BACKGROUND_GRADIENT=False
WITH_TOPOGRAPHY=False
WITH_RIVERS_AND_BOUNDARIES=True
WITH_DATETIME=True

RESUME_MODE=bool(RESUME_AT_FILE!="")
if RESUME_MODE:
    RESUME_AT_FILE = os.path.split(RESUME_AT_FILE)[1]
    print "Resuming at file " + RESUME_AT_FILE

# Color tables

# IR 10.8
#msevi_l15_ir_108_ctable = "GnBu"
msevi_l15_ir_108_ctable = "xray"
invert_msevi_l15_ir_108_ctable = 1

# C-Band RX
#cband_radolan_rx_ctable = "RdYlBu"
#invert_cband_radolan_rx_ctable = 1
#cband_radolan_rx_ctable = "Spectral"
cband_radolan_rx_ctable = "hot_desaturated"
invert_cband_radolan_rx_ctable = 0

# LINET
linet_oase_tl_ctable = "contoured"
invert_linet_oase_tl_ctable = 0

# Conversion program params
CONVERSION_PARAMS  = "-t cluster "
CONVERSION_PARAMS += " -v msevi_l15_ir_108"
CONVERSION_PARAMS += " --write-as-xml=false"
CONVERSION_PARAMS += " --extract-skin=false"
CONVERSION_PARAMS += " --vtk-dimensions x,y"

#-f herz-oase-20120805t1510utc-0500m-bonnjue-3d-v01a_clusters.nc "

# binaries
DYLD_LIBRARY_PATH="/usr/local/lib"
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
conversion_bin = bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk"
print "Conversion Command:"
print conversion_bin

# Silent
SuppressMessages(True)
SuppressQueryOutputOn()

# Set view and annotation attributes

a = GetAnnotationAttributes()
a.databaseInfoFlag = 0
a.timeInfoFlag = 0
a.userInfoFlag = 0
a.legendInfoFlag = 1
a.axes2D.visible=1
a.axes2D.xAxis.title.visible=0
a.axes2D.yAxis.title.visible=0
a.axes2D.autoSetScaling=0

if OVERLAY:
    a.legendInfoFlag=0
    a.axes2D.visible=0

SetAnnotationAttributes(a)

print ColorTableNames()

# Add gray/black background gradient
if WITH_BACKGROUND_GRADIENT:
    visitUtils.add_background_gradient();

if not RESUME_MODE:
    print "Cleaning up *.vtk *.vtr *.png"
    return_code=call("rm -f *.vtk *.vtr *.png", shell=True)
else:
    print "Cleaning up *.vtk *.vtr"
    return_code=call("rm -f *.vtk *.vtr", shell=True)

# Set view to nationwide composite
visit2D.set_view_to_radolan();

print "-- Creating colortables ---"
num_colors = visitUtils.create_cluster_colortable("cluster_colors")

if WITH_TOPOGRAPHY:
    visitUtils.create_topography_colortable()

print "    done."

# Glob the netcdf directory
netcdf_files = glob.glob(NETCDF_DIR+"/*.nc");

for netcdf_file in netcdf_files:

    # construct the cluster filename and find it
    # in the cluster directory
    
    netcdf_path,filename    = os.path.split(netcdf_file);
    basename                = os.path.splitext(filename)[0]

    cluster_file            = CLUSTER_DIR+"/"+basename+"-clusters.nc"
    label_file              = basename+"-clusters_centers.vtk"

    #print "netcdf_file  = " + netcdf_file
    #print "filename     = " + filename
    #print "basename     = " + basename
    #print "cluster_file = " + cluster_file
    #print "label_file   = " + label_file
    
    # If we're resuming, skip until the appointed file
    
    if RESUME_MODE:
        if filename == RESUME_AT_FILE:
            RESUME_MODE=False
        else:
            continue
    
    # check if the files both exist
    print "Processing file "+netcdf_file+" and cluster file "+cluster_file
    if not os.path.exists(cluster_file):
        print "Cluster file does not exist. Skipping."
        continue
    
    if (CREATE_SOURCE_MOVIE or (CREATE_CLUSTERS_MOVIE and not OVERLAY_CLUSTERS)):

        if not OVERLAY:

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

        # now plot the data
        OpenDatabase(netcdf_file);

        num_base_plots = GetNumPlots();

        # IR 10.8

        if not OVERLAY:

            visit2D.add_pseudocolor(netcdf_file,"msevi_l15_ir_108",msevi_l15_ir_108_ctable,1,1)
            p = PseudocolorAttributes()
            p.colorTableName=msevi_l15_ir_108_ctable
            p.minFlag,p.maxFlag=1,1
            p.min,p.max=10,71.9493
            p.invertColorTable = invert_msevi_l15_ir_108_ctable
            SetActivePlots(num_base_plots+1)
            SetPlotOptions(p)
            
            AddOperator("Threshold")
            t = ThresholdAttributes();
            t.upperBounds=(71.9493)
            SetOperatorOptions(t)

        # C-band RX

        visit2D.add_pseudocolor(netcdf_file,"cband_radolan_rx",cband_radolan_rx_ctable,1,1)
        p = PseudocolorAttributes()
        p.colorTableName=cband_radolan_rx_ctable
        p.minFlag,p.maxFlag=1,1
        p.min,p.max=0,65
        p.invertColorTable=invert_cband_radolan_rx_ctable
        SetActivePlots(num_base_plots+2)
        SetPlotOptions(p)
        AddOperator("Threshold")
        t = ThresholdAttributes();
        t.lowerBounds=(15)
        SetOperatorOptions(t)

        if OVERLAY:
            DrawPlots()
            visitUtils.save_window("cband_rx_",1)
            DeleteAllPlots()

        # LINET

        visit2D.add_pseudocolor(netcdf_file,"linet_oase_tl",linet_oase_tl_ctable,1,1)
        p = PseudocolorAttributes()
        p.colorTableName=linet_oase_tl_ctable
        p.minFlag,p.maxFlag=1,1
        p.min,p.max=1,20
        p.invertColorTable = invert_linet_oase_tl_ctable
        SetActivePlots(num_base_plots+3)
        SetPlotOptions(p)
        AddOperator("Threshold")
        t = ThresholdAttributes();
        t.lowerBounds=(1)
        SetOperatorOptions(t)

        DrawPlots()
            
        if OVERLAY:
            visitUtils.save_window("linet_",1)
            DeleteAllPlots()

        if not OVERLAY:
            visitUtils.save_window("source_",1)

    if CREATE_CLUSTERS_MOVIE:
        
        if not OVERLAY_CLUSTERS:
            DeleteAllPlots()
            ClearWindow()

        start_time = time.time()
        print "-- Converting clusters to .vtr --"

        # build the clustering command
        command=conversion_bin+" -f "+cluster_file+" "+CONVERSION_PARAMS
        print command
        return_code = call( command, shell=True)

        print "    done. (%.2f seconds)" % (time.time()-start_time)
        print "-- Rendering cluster scene --"
        start_time = time.time()

        if not (OVERLAY_CLUSTERS and OVERLAY):
    
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

        # Add the clusters
        basename = CLUSTER_DIR+"/"
        visit2D.add_clusters_with_colortable(basename,"_cluster_","cluster_colors",num_colors)

        # Add modes as labels
        visitUtils.add_labels(label_file,"geometrical_center")

        # save as image
        DrawPlots()

        if OVERLAY:
            visitUtils.save_window("clusters_",1)
        else:
            visitUtils.save_window("tracking_",1)

        DeleteAllPlots()
        ClearWindow()

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


if OVERLAY:

    # run imagemagick
    background=M3D_HOME+"/visit/1024x1024-transparent.png"

    # cband

    files=sorted(glob.glob("cband_rx_*.png"))
    for f in files:
        return_code=call("convert "+f+" -transparent white -crop 765x765+207+104 "+f, shell=True)
        return_code=call("convert -composite -geometry +205+105 "+background+" "+f+" "+f, shell=True)

    # linet

    files=sorted(glob.glob("linet_*.png"))
    for f in files:
        return_code=call("convert "+f+" -transparent white -crop 765x765+207+104 "+f, shell=True)
        return_code=call("convert -composite -geometry +205+105 "+background+" "+f+" "+f, shell=True)

    # clusters

    files=sorted(glob.glob("clusters_*.png"))
    for f in files:
        return_code=call("convert "+f+" -transparent white -crop 765x765+207+104 "+f, shell=True)
        return_code=call("convert -composite -geometry +205+105 "+background+" "+f+" "+f, shell=True)

    # create loops
    if CREATE_SOURCE_MOVIE:
        
        visitUtils.create_movie("cband","cband_rx.gif")
        visitUtils.create_movie("cband","cband_rx.m4v")

        visitUtils.create_movie("linet","linet.gif")
        visitUtils.create_movie("linet","linet.m4v")

    if  CREATE_CLUSTERS_MOVIE:
        
        visitUtils.create_movie("clusters","clusters.gif")
        visitUtils.create_movie("clusters","clusters.m4v")

else:

    # create loops
    if CREATE_SOURCE_MOVIE:
        visitUtils.create_movie("source_","source.gif")
        visitUtils.create_movie("source_","source.m4v")

    if  CREATE_CLUSTERS_MOVIE:
        visitUtils.create_movie("tracking_","tracking.gif")
        visitUtils.create_movie("tracking_","tracking.m4v")

# Clean up
print "Cleaning up ..."
return_code=call("mkdir images", shell=True)
return_code=call("mv *.png images", shell=True)
return_code=call("mkdir movies", shell=True)
return_code=call("mv *.gif *.m4v movies", shell=True)
return_code=call("rm -f *.vt* visitlog.py", shell=True)

quit()
