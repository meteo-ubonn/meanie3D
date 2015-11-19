'''
The MIT License (MIT)

(c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import glob
import os
import subprocess
import sys
import time
import visit

# Own packages
import meanie3D.visualisation
import utils

from meanie3D.app import external
external.locateCommands(['meanie3D-cfm2vtk'])

# ------------------------------------------------------------------------------
# Adds clusters with names "*_infix_*.vtk" to the current visualisation window.
# The clusters are colored according to the cluster_colors colortable,
# which is created at the beginning of the visualisation code.
#
# \param infix (e.g. "_untracked_clusters_")
# \param color_table_name Name of the cluster color table
# \param color_count Number of colors in the cluster color table
# \param configuration Options to overwrite the options the clusters are
# plotted with. Clusters are added as pseudocolor plots. For information
# on the available options, check PseudocolorAttributes()
#
# ------------------------------------------------------------------------------
def addClusters(infix,configuration):
    # Get the configuration
    clusterOptions = utils.getValueForKeyPath(configuration,'postprocessing.clusters.visit.cluster')
    if not clusterOptions:
        sys.stderr.write("ERROR:no value found for postprocessing.clusters.visit.cluster")
        return

    # now the clusters
    cluster_pattern = "*"+infix+"*.vt*"
    print "Looking for cluster files at " + cluster_pattern
    cluster_list = sorted(glob.glob(cluster_pattern))
    if not cluster_list:
        # all?
        cluster_list = sorted(glob.glob("*_clusters_all.vt*"))
    print "Found %d cluster files." % len(cluster_list)
    for cluster_file in cluster_list:
        utils.addPseudolorPlot(cluster_file,clusterOptions)

    return

# ------------------------------------------------------------------------------
# Generic routine for visualizing 3D clusters in two perspectives
# \param: configuration
# ------------------------------------------------------------------------------
def run(conf):
    #pp = pprint.PrettyPrinter()
    #pp.pprint(conf)

    # Make sure the global configuration is in place
    utils.runGlobalVisitConf(conf)
    #pdb.set_trace()

    clusterConf = utils.getValueForKeyPath(conf,'postprocessing.clusters')
    if not clusterConf:
        print "No configuration for cluster postprocessing. Nothing to do."
        return 0

    visitConf = utils.getValueForKeyPath(clusterConf,'visit')
    if not visitConf:
        print "No configuration for visuals. Nothing to do."
        return 0

    views = utils.getValueForKeyPath(visitConf,'views')

    # Set up background gradient, axis labels etc.
    utils.setAnnotations(visitConf,'annotationAttributes')

    # Set view and annotation attributes
    utils.setAnnotations(conf,'postprocessing.clusters.visit.annotations')

    if not utils.getValueForKeyPath(conf,'resume'):
        print "Removing results from previous runs"
        subprocess.call("rm -rf images movies *.vtk *.vtr *tracking*.png *source*.png", shell=True)
    else:
        print "Removing intermediary files from previous runs"
        subprocess.call("rm -f *.vtk *.vtr", shell=True)

    # Glob the netcdf directory
    print "Current work directory: " + os.path.abspath(os.getcwd())
    print "Processing files in directory " + conf['source_directory']
    netcdf_files = sorted(glob.glob(conf['source_directory']+"/*.nc"))

    # Keep track of number of images to allow
    # forced re-set in time to circumvent the
    # Visit memory leak
    image_count=0


    for netcdf_file in netcdf_files:
        # construct the cluster filename and find it
        # in the cluster directory
        netcdf_path,filename    = os.path.split(netcdf_file);
        basename                = os.path.splitext(filename)[0]
        cluster_file            = conf['cluster_directory']+"/"+basename+"-clusters.nc"
        label_file              = basename+"-clusters_centers.vtk"
        displacements_file      = basename+"-clusters_displacements.vtk"

        print "netcdf_file  = " + netcdf_file
        print "cluster_file = " + cluster_file

        # check if the files both exist
        if not os.path.exists(cluster_file):
            print "Cluster file does not exist. Skipping."
            continue

        # predict the filenames for checking on resume
        number_postfix = str(image_count).rjust(4,'0')

        source_open = False
        skip_source = False

        if conf['resume'] == True:
            exists = utils.images_exist(visitConf['views'],"source",image_count)
            if exists == "all":
                print "Source visualization "+number_postfix+" exists. Skipping."
                skip_source = True
            elif exists == "partial":
                print "Deleting partial visualization " + number_postfix
                utils.delete_images(conf,"source",image_count)

        if skip_source == False:

            if utils.getValueForKeyPath(clusterConf,'createSourceMovie'):

                # Add ancillary background data
                utils.plotMapdata(visitConf,'map')

                # Add timestamp
                if utils.getValueForKeyPath(clusterConf,'showDateTime'):
                    utils.add_datetime(netcdf_file)

                # Add source data and threshold it
                print "Plotting source data ..."
                start_time = time.time()

                utils.addPseudocolorPlots(netcdf_file,visitConf,'source.plots')
                source_open = True
                visit.DrawPlots()

                utils.saveImagesForViews(views,"source")

                visit.DeleteAllPlots()
                visit.ClearWindow()

                print "    done. (%.2f seconds)" % (time.time()-start_time)

        if utils.getValueForKeyPath(clusterConf,'visualiseClusters'):

            skip = False
            if conf['resume'] == True:

                exists = utils.images_exist(visitConf['views'],"tracking",image_count)
                if exists == "all":
                    print "Cluster visualization "+number_postfix+" exists. Skipping."
                    skip = True
                elif exists == "partial":
                    print "Deleting partial cluster visualization " + number_postfix
                    utils.delete_images(conf,"tracking",image_count)

            if skip == False:

                # Run the conversion
                print "-- Converting clusters to .vtr --"
                start_time = time.time()
                params = "-f %s %s" \
                         % (cluster_file,utils.getValueForKeyPath(conf,'postprocessing.clusters.meanie3D-cfm2vtk'))
                if utils.getValueForKeyPath(conf,'postprocessing.clusters.showDisplacementVectors'):
                    params.append(" --write-displacement-vectors")

                if utils.getValueForKeyPath(conf,'data.vtkDimensions'):
                    params.append(" --vtk-dimensions=%s" % conf['vtkDimensions'])

                # pdb.set_trace();
                meanie3D.app.external.execute_command('meanie3D-cfm2vtk', params)
                print "    done. (%.2f seconds)" % (time.time()-start_time)

                print "-- Rendering cluster scene --"
                start_time = time.time()

                # Add ancillary background data
                utils.plotMapdata(visitConf,'map')

                # Add timestamp
                if utils.getValueForKeyPath(clusterConf,'showDateTime'):
                    utils.add_datetime(netcdf_file)

                # Add background source data
                if utils.getValueForKeyPath(clusterConf,'showSourceBackground'):
                    utils.addPseudocolorPlots(netcdf_file,visitConf,'sourceBackground.plots')
                    source_open = True

                # Add the clusters
                basename = conf['cluster_directory'] + "/"
                addClusters("_cluster_",conf)

                # Add modes as labels
                labelConf = utils.getValueForKeyPath(visitConf,'label')
                if labelConf:
                    utils.addLabelPlot(label_file,labelConf)

                # Add displacement vectors
                if utils.getValueForKeyPath(clusterConf,'showDisplacementVectors'):
                    vectorConf = utils.getValueForKeyPath(visitConf,'displacementVectors')
                    if vectorConf:
                        utils.addVectorPlot(displacements_file,vectorConf)

                visit.DrawPlots()
                utils.saveImagesForViews(views,"tracking")

                print "    done. (%.2f seconds)" % (time.time()-start_time)

        # clean up

        visit.DeleteAllPlots();
        visit.ClearWindow()
        if source_open:
            visit.CloseDatabase(netcdf_file)
            visit.CloseDatabase(label_file)
        utils.close_pattern(basename+"*.vtr")
        utils.close_pattern(basename+"*.vtk")
        subprocess.call("rm -f *.vt*", shell=True)

        # periodically kill computing engine to
        # work around the memory leak fix
        image_count=image_count+1;

        # TODO: check if this is still necessary and re-implement if it is.
        #if image_count % 100 == 0:
        #    visit.CloseComputeEngine()

    # close mapstuff
    # TODO: might need to keep track of open database files and
    # close them.

    # Use imagemagick to use image sequences to make movies

    movieFormats = utils.getValueForKeyPath(clusterConf,'movieFormats')
    if utils.getValueForKeyPath(clusterConf,'createSourceMovie'):
        utils.createMoviesForViews(views,"source", movieFormats)

    if utils.getValueForKeyPath(clusterConf,'createClusterMovie'):
        utils.createMoviesForViews(views,"tracking", movieFormats)

    # clean up
    print "Cleaning up ..."
    subprocess.call("mkdir images", shell=True)
    subprocess.call("mv *tracking_*.png images", shell=True)
    subprocess.call("mv *source_*.png images", shell=True)
    subprocess.call("mkdir movies", shell=True)
    subprocess.call("mv *source*.gif *tracking*.gif *.m4v movies", shell=True)
    subprocess.call("rm -f *.vt* visitlog.py", shell=True)
    return

