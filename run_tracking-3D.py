3#!/usr/bin/python
import glob
import sys
import os
from subprocess import call

# python script template for tracking
# TODO: clustering parameters should be passed in from the bash script

MEANIE3D_HOME     = "M3D_HOME"
DYLD_LIBRARY_PATH = "DL_PATH"
NETCDF_DIR        = "SOURCE_DIR"
CLUSTERING_PARAMS = "-d z,y,x -v zh -w zh -r 3,7,7,100 --drf-threshold 0.75 -s 12 -t 10 --write-variables-as-vtk=zh --vtk-dimensions=x,y,z" 
TRACKING_PARAMS   = "-t zh --verbosity 3 --write-vtk --vtk-dimensions=x,y,z --wr=1.0 --ws=1.0 --wt=0.0"

# print parameters

print "Running clustering on directory "+NETCDF_DIR
print "MEANIE3D_HOME="+MEANIE3D_HOME
print "DYLD_LIBRARY_PATH="+DYLD_LIBRARY_PATH

toggled_maintain=False

# delete any previous results
# results are stored in the work directory

# delete previous results
return_code=call("rm *-clusters.nc", shell=True)
return_code=call("rm *_clusters_*.vtk", shell=True)

# binaries
bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
detection_bin = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-detect"
tracking_bin  = bin_prefix + "M3D_HOME" + "/Debug/meanie3D-track"

print "Detection Command:"
print detection_bin
print "Tracking Command:"
print tracking_bin

# Cluster color tables
col_tables = ["Purples","Blues","Oranges","Greens","Reds","Set1","Set2","Set3"]

# Set view and annotation attributes
a = GetAnnotationAttributes()
a.axes3D.visible=1
a.axes3D.setBBoxLocation=1
a.axes3D.bboxLocation = (-366, -111, -4340, -4105, 0, 15)
a.axes3D.autoSetScaling=0
a.userInfoFlag=0
a.timeInfoFlag=0
a.legendInfoFlag=0
a.databaseInfoFlag=1
SetAnnotationAttributes(a)
print a

# Modify view parameters
v = GetView3D()
v.focus=(-238.5,-4222.5,50.0)
SetView3D(v)

s = GetSaveWindowAttributes()
print s

# Get a list of the files we need to process
netcdf_pattern = NETCDF_DIR + "/*.nc"
netcdf_list=glob.glob(netcdf_pattern)
last_cluster_file=""
file_count=0

# Process the files one by one
for netcdf_file in netcdf_list:

    print "Processing " + netcdf_file

    basename = os.path.basename(netcdf_file)
    cluster_file=os.path.splitext(basename)[0]+"-clusters.nc"
    vtk_file=os.path.splitext(basename)[0]+".vtk"

    # build the clustering command
    command=detection_bin+" -f "+netcdf_file+" -o "+cluster_file + " " + CLUSTERING_PARAMS
    
    # if it's the first file, write the clusters out from here
    # else from the tracking tool
    if last_cluster_file == "" :
        command = command + " --write-clusters-as-vtk"
        
    # execute
    print command
    return_code = call( command, shell=True)

    # if we have a previous scan, run the tracking command

    if last_cluster_file != "" :
        command=tracking_bin+" -p "+last_cluster_file+" -c "+cluster_file+" " + TRACKING_PARAMS
        print command
        return_code = call( command, shell=True)

    # keep track
    last_cluster_file=cluster_file
 
    # open the file and add the plot
    OpenDatabase(vtk_file)
    AddPlot("Pseudocolor", "zh")

    p = PseudocolorAttributes()
    p.colorTableName = "hot_desaturated"
    #p.legendFlag=1
    #p.lightingFlag=0
    #p.invertColorTable=0
    #p.minFlag,p.maxFlag = 1,1
    #p.min,p.max = 0.0, 50.0
    #p.opacity=1.0
    SetPlotOptions(p)

    # save the source file
    s = GetSaveWindowAttributes()
    s.fileName="source_";
    s.progressive=1
    print s
    SetSaveWindowAttributes(s)
    
    # Draw, save and remove
    DrawPlots()
    SaveWindow()
    DeleteAllPlots();

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
        cp.opacity=1
        SetPlotOptions(cp)

    # Get it all plotted
    DrawPlots()

    if toggled_maintain != True :
        ToggleMaintainViewMode()
        toggled_maintain=True

    # save the source file

    s = GetSaveWindowAttributes()
    s.fileName="tracking_";
    s.progressive=1
    print s
    SetSaveWindowAttributes(s)

    DrawPlots()
    SaveWindow()
    DeleteAllPlots();

    # clean up

    ClearWindow()
    CloseDatabase(vtk_file)
    for cluster_file in cluster_list:
        CloseDatabase(cluster_file)
    #return_code=call("rm *.vtk", shell=True)

    file_count += 1

print "Done. Closing Visit."
Close();
