#!/usr/bin/python
import glob
import sys
import os
from subprocess import call

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
DYLD_LIBRARY_PATH = "DL_PATH"
NETCDF_FILE       = "SOURCE_FILE"

# Derive some useful filenames

basename = os.path.basename(NETCDF_FILE)
vtk_file=os.path.splitext(basename)[0]+".vtk"

# Command parameters

# TODO: infer the correct set from the command line
# parameter automatically

# OASE Komposit
#CLUSTERING_PARAMS = "-d z,y,x -v zh -w zh -r 5,10,10,100 --drf-threshold 0.75 -s 16 -t 10 --write-variables-as-vtk=zh --vtk-dimensions=x,y,z" 

# RADOLAN
CLUSTERING_PARAMS = "-v reflectivity -w reflectivity -d x,y -r 5,5,200 -t 10 -m 10 --write-clusters-as-vtk"

# print parameters

print "Running clustering on file "+NETCDF_FILE
print "MEANIE3D_HOME="+MEANIE3D_HOME
print "DYLD_LIBRARY_PATH="+DYLD_LIBRARY_PATH

toggled_maintain=False

# delete any previous results
# results are stored in the work directory
return_code=call("rm *cluster_*.vtk", shell=True)

# binaries
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
detection_bin = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-detect"
tracking_bin  = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-track"

print "Detection Command:"
print detection_bin
print "Tracking Command:"
print tracking_bin

# scales to use
#scales = ["4","8","16","32","64","128","256","512","1024"]
scales = ["10","25","50","75"]
resume_at_scale=""

# DRFs to use
drfs = ["0.0","0.5","0.75","1.0"]
resume_at_drf=""

# Cluster color tables
col_tables = ["Purples","Blues","Oranges","Greens","Reds"]

# Suppress dialogues
SuppressQueryOutputOn()

# Set view and annotation attributes
a = GetAnnotationAttributes()
a.axes2D.visible=1
a.axes2D.autoSetScaling=0
a.userInfoFlag=0
a.timeInfoFlag=0
a.legendInfoFlag=0
a.databaseInfoFlag=1
SetAnnotationAttributes(a)

# Modify view parameters
v = GetView3D()
v.focus=(-238.5,-4222.5,50.0)
SetView3D(v)

# flag for avoiding creating the source image and 
# source vtk files multiple times
source_saved=False


for scale in scales:

    if resume_at_scale!="" and resume_at_scale!=scale:
        continue
    else:
        resume_at_scale=""

    for drf in drfs:

        if resume_at_drf!="" and resume_at_drf!=drf:
            continue
        else:
            resume_at_drf=""

        PARAMS = CLUSTERING_PARAMS + " --drf-threshold " +drf + " -s " + scale

        if source_saved == False:
           #radolan case
           PARAMS = PARAMS + " --write-variables-as-vtk=reflectivity" 

        # build the clustering command
        cluster_file=os.path.splitext(basename)[0]+"_clusters_scale_"+scale+"_drf_"+drf+".nc"
        command=detection_bin+" -f "+NETCDF_FILE+" -o "+cluster_file + " " + PARAMS
    
        # execute
        print command
        return_code = call( command, shell=True)

        # open the file and add the plot
        OpenDatabase(vtk_file)
        AddPlot("Pseudocolor", "reflectivity")

        p = PseudocolorAttributes()
        p.colorTableName = "xray"
        #p.legendFlag=1
        p.lightingFlag=1
        #p.invertColorTable=1
        p.minFlag,p.maxFlag = 1,1
        p.min,p.max = 0.0, 50.0
        p.opacity=1.0
        SetPlotOptions(p)

        if source_saved==False:
            # save the source file
            s = GetSaveWindowAttributes()
            s.fileName="source_"
            s.progressive=0
            SetSaveWindowAttributes(s)
            DrawPlots()
            SaveWindow()
            DeleteAllPlots();
            source_saved=True

        # now the clusters
        cluster_pattern=basename+"*cluster*.vtk"
        print "Looking for cluster files at " + cluster_pattern
        cluster_list=glob.glob(cluster_pattern)

        print "Processing clusters:"
        print cluster_list

        for cluster_file in cluster_list:
        
            # figure out cluster number for choice of color table
            # extract the number from the cluster filename
            # to select the color table based on the ID
        
            cluster_num=0

            try:
                print "Extracting cluster number from " + cluster_file;
                start = cluster_file.index("_cluster_",0) + len("_cluster_")
                end = cluster_file.index(".vtk",start)
                cluster_num = int( cluster_file[start:end] );
                print "Plotting cluster #" + str(cluster_num)
            except ValueError:
                print "Illegal filename " + fname
                continue

            OpenDatabase(cluster_file)
            AddPlot("Pseudocolor", "cluster")
            cp=PseudocolorAttributes();
            #cp.minFlag,cp.maxFlag = 1,1
            #cp.min,cp.max = 10.0, 60.0
            cp.pointSizePixels=5
            cp.legendFlag=0 
            cp.lightingFlag=1
            cp.invertColorTable=0
            index = cluster_num % len(col_tables);
            cp.colorTableName = col_tables[index];
            cp.opacity=0.05
            SetPlotOptions(cp)

        # Get it all plotted
        DrawPlots()

        if toggled_maintain != True :
            ToggleMaintainViewMode()
            toggled_maintain=True

        # Save the image
        s = GetSaveWindowAttributes()
        s.fileName=basename+"_clusters_scale_"+scale+"_drf_"+drf+"_"
        s.progressive=0
        SetSaveWindowAttributes(s)
        DrawPlots()
        SaveWindow()
        DeleteAllPlots();

        # clean up
        ClearWindow()
        CloseDatabase(vtk_file)
        for cluster_file in cluster_list:
            CloseDatabase(cluster_file)
        return_code=call("rm *cluster_*.vtk", shell=True)

print "Done. Closing Visit."
Close();
