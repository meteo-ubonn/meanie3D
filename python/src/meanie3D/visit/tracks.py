#!/usr/bin/python
__author__ = "juergen.simon@uni-bonn.de"

# ------------------------------------------------------------------------------
# Filename: tracks_visit.py
#
# This module bundles python routines for running inside Visit.
#
# @author Juergen Simon (juergen.simon@uni-bonn.de)
# ------------------------------------------------------------------------------

import meanie3D.app.utils

##
# Sets the view parameters on a 2D or 3D view. The conf dictionary
# may contain all keys that are used on an object  returned by
# visit's GetView2D() or GetView3D() call.
#
# \param conf
# \param dimensionality (2 or 3)
#
def setView(conf,dimensions):
    if (dimensions == 2):
        v = GetView2D();
        v.windowCoords = tuple(conf['windowCoords'])
        v.viewportCoords = tuple(conf['viewportCoords'])
        SetView2D(v)
    else:
        v = GetView3D();
        v.viewNormal = tuple(conf["viewNormal"])
        v.focus = tuple(conf["focus"])
        v.viewUp = tuple(conf["viewUp"])
        v.viewAngle = conf["viewAngle"]
        v.parallelScale = conf["parallelScale"]
        v.nearPlane = conf["nearPlane"]
        v.farPlane = conf["farPlane"]
        v.imagePan = tuple(conf["imagePan"])
        v.imageZoom = conf["imageZoom"]
        v.perspective = conf["perspective"]
        v.eyeAngle = conf["eyeAngle"]
        v.centerOfRotationSet = conf["centerOfRotationSet"]
        v.centerOfRotation = tuple(conf["centerOfRotation"])
        v.axis3DScaleFlag = conf["axis3DScaleFlag"]
        v.axis3DScales = tuple(conf["axis3DScales"])
        v.shear = tuple(conf["shear"])
        SetView3D(v)
##
# Plots mapdata according to the configuration given. The configuration
# has the follwing format (in json):
#
# "map" : {
#   "mapDataFile" : "$MEANIE3D_HOME/oase-mapdata/oase-mapdata.nc",
#   "variables" : ["as_zonal/national_boundaries_2D","as_zonal/national_rivers_2D"],
#   "min" : [0,0],
#   "max" : [1,1],
#   "colorTableName" : ["Greys","hot"],
#   "invertColorTable" : [0,1],
#   "opacities" : [1,1]
# }
#
# Note that $ variables will be replaced with the value found in the
# environment.
#
def plotMapdata(conf):
    if (conf['mapDataFile']):
        # Expand dollar variables
        mapFile = os.path.expandvars(conf['mapDataFile'])
        # Check if data file exists
        if not os.path.exists(mapFile):
            print "ERROR:could not find map file at " + mapFile
            return

        # open the file and add the plot
        OpenDatabase(mapFile)

        # Iterate over variables
        for i in range(0,len(conf['variables'])):
            variable = conf['variables'][i]
            AddPlot("Pseudocolor", variable)
            p = PseudocolorAttributes()
            p.lightingFlag = 0
            p.legendFlag = 0;
            p.colorTableName = conf['colorTableName'][i]
            p.invertColorTable = conf['invertColorTable'][i];
            p.opacity = conf['opacity'][i]
            p.minFlag,p.maxFlag = 1,1
            p.min,p.max = conf['min'][i], conf['max'][i]
            SetPlotOptions(p)

        DrawPlots()

##
# Plots a list of tracks from .vtk files produced by
# meanie3D-trackstats --write-center-tracks-as-vtk.
# \param conf configuration dictionary containing the following keys
#
# 'cluster_directory' : directory with the cluster results
# 'meanie3d_home' : home directory of meanie3D (for the mapstuff file and modules)
#
#  "visit" : {
#       "dimensions" : 2,
#       "view" : {
#           ...
#       },
#       "showBackgroundGradient" : false,
#       "map" : {
#           ...
#       },
#       "tracks" : {
#           ...
#       }
#  }
#
def plotTracks(conf):

    visitConf = meanie3D.app.utils.getValueForKeyPath(conf,'tracks.visit')
    if not visitConf:
        print "No configuration for visuals. Nothing to do."
        return 0

    # Silent
    SuppressMessages(2)
    SuppressQueryOutputOn()

    # Check the configuration
    if (not visitConf['dimensions']):
        print "ERROR:configuration must contain 'dimensions'"
        return -1
    if (not visitConf['dimensions'] in [2,3]):
        print "ERROR:'dimensions' can only be 2 or 3."
        return -1

    # Set annotation attributes.
    # TODO: find a sensible method of exposing this to config
    a = GetAnnotationAttributes()
    if (visitConf['dimensions'] == 2):
        # 2D
        a.axes2D.visible=1
        a.axes2D.autoSetScaling=0
        a.axes2D.xAxis.title.visible=0
        a.axes2D.yAxis.title.visible=0
    else:
        # 3D
        a.axes3D.visible=1
        a.axes3D.autoSetScaling=0
        a.axes3D.xAxis.title.visible=0
        a.axes3D.yAxis.title.visible=0
        a.axes3D.zAxis.title.visible=0

    a.legendInfoFlag=1
    a.databaseInfoFlag=0
    a.userInfoFlag=0
    a.timeInfoFlag=0
    SetAnnotationAttributes(a)

    # Add gray/black background gradient
    #if (visitConf['showBackgroundGradient'] == True):
        # TODO: show background gradient
        # visitUtils.add_background_gradient();

    # Set view to nationwide composite
    if (visitConf['view']):
        setView(visitConf['view'],visitConf['dimensions'])

    if (visitConf['map']):
        plotMapdata(visitConf['map'])

    if (visitConf['track']):
        trackConf = visitConf['track']

        # Plot the Tracks
        # track_pattern = conf['tracks_dir'] + "/*-track_*.vtk"
        os.chdir(conf['tracks_dir'])
        track_pattern = "*-track_*.vtk"
        list = sorted(glob.glob(track_pattern))

        print "Looking with pattern " + track_pattern
        print "Found %d track files." % len(list)
        count = 0;
        for fname in list:

            # add plot
            OpenDatabase(conf['tracks_dir'] + "/" + fname);
            AddPlot("Pseudocolor",trackConf['variable'])

            # Configured attributes
            cp=PseudocolorAttributes();
            cp.pointSizePixels = trackConf['pointSizePixels']
            cp.minFlag,cp.maxFlag = trackConf['minFlag'],trackConf['maxFlag']
            cp.min,cp.max = trackConf['min'],trackConf['max']
            cp.opacity = trackConf['opacity']
            cp.colorTableName = trackConf['colorTableName']
            cp.invertColorTable = trackConf['invertColorTable']

            # plot the legend for the first one
            if count==0:
                cp.legendFlag=trackConf['legendFlag']
            else:
                cp.legendFlag=0
            SetPlotOptions(cp)
            count = count+1

            # in case the script is being debugged, exit the script
            # after 10 tracks. This could be configured
            if conf['tracks']['debugVisitScript'] and count > 10:
                quit()

        DrawPlots()
        utils.save_window("tracks",0)

    quit()
