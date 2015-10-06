#!/usr/bin/python
__author__ = "juergen.simon@uni-bonn.de"

# ------------------------------------------------------------------------------
# Filename: tracks_visit.py
#
# This module bundles python routines for running inside Visit.
#
# @author Juergen Simon (juergen.simon@uni-bonn.de)
# ------------------------------------------------------------------------------

# System packages
import os
import pprint

# Visit package
from visit import *

# Own packages
from .utils import *
from meanie3D.app.utils import getValueForKeyPath
from meanie3D.app.utils import getSafe


##
# Plots a list of tracks from .vtk files produced by
# meanie3D-trackstats --write-center-tracks-as-vtk.
#
# \param:conf Configuration dictionary
#
def plotTracks(conf):

    visitConf = getValueForKeyPath(conf,'postprocessing.tracks.visit')
    if not visitConf:
        print "No configuration for visuals. Nothing to do."
        return 0

    # Silent
    # SuppressMessages(2)
    # SuppressQueryOutputOn()

    # Check the configuration
    dimensions = getValueForKeyPath(conf,'data.dimensions')

    if (not dimensions):
        print "ERROR:configuration must contain 'dimensions'"
        return -1
    if (not len(dimensions) in [2,3]):
        print "ERROR:Can only process 2D or 3D"
        return -1

    # Set up background gradient, axis labels etc.
    setAnnotations(conf,'postprocessing.tracks.visit.annotationAttributes')

    # Set the view straight
    setView(conf,'postprocessing.tracks.visit.view')

    # Plot the map data
    plotMapdata(conf,'postprocessing.tracks.visit.map')

    # Plot the tracks
    trackPlotConf = getValueForKeyPath(conf,'postprocessing.tracks.visit.track')
    if trackPlotConf:

        # Save value of legend flag
        legendFlag = getValueForKeyPath(trackPlotConf,'PseudocolorAttributes.legendFlag')

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
            trackFile = conf['tracks_dir'] + "/" + fname

            # plot the legend for the first one only
            if (count == 1) and legendFlag:
                meanie3D.app.utils.setValueForKeyPath(trackPlotConf,'PseudocolorAttributes.legendFlag',0)

            # Plot the actual track data
            add_pseudocolor(trackFile,trackPlotConf)

            count = count + 1

            # in case the script is being debugged, exit the script
            # after 10 tracks. This could be configured
            # if getValueForKeyPath(conf,'postprocessing.debugVisitScript') and count > 10

        # Restore flag value
        meanie3D.app.utils.setValueForKeyPath(trackPlotConf,'PseudocolorAttributes.legendFlag', legendFlag)

    DrawPlots()
    save_window("tracks",0)
    return
