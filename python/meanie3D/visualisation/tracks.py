##
# This modules contains code to visualise tracks in visit.
#

import glob
import os
import pprint
import visit
import utils

##
# Plots a list of tracks from .vtk files produced by
# meanie3D-trackstats --write-center-tracks-as-vtk.
#
# \param:conf Configuration dictionary
#
def run(conf):

    pp = pprint.PrettyPrinter()
    # pp.pprint(conf)

    # Make sure the global configuration is in place
    utils.runGlobalVisitConf(conf)

    visitConf = utils.getValueForKeyPath(conf,'postprocessing.tracks.visit')
    if not visitConf:
        print "No configuration for visuals. Nothing to do."
        return 0

    # Set up background gradient, axis labels etc.
    utils.setAnnotations(conf,'postprocessing.tracks.visit.annotationAttributes')

    # Set the view straight
    utils.setView(conf,'postprocessing.tracks.visit.view')

    # Plot the map data
    utils.plotMapdata(conf,'postprocessing.tracks.visit.map')

    # Plot the tracks
    trackPlotConf = utils.getValueForKeyPath(conf,'postprocessing.tracks.visit.track')
    # pp.pprint(trackPlotConf)

    currentDirectory = os.path.abspath(os.getcwd())
    os.chdir(conf['tracks_dir'])

    if trackPlotConf:

        # Save value of legend flag
        legendFlag = trackPlotConf['PseudocolorAttributes']['legendFlag']

        # Plot the Tracks
        # track_pattern = conf['tracks_dir'] + "/*-track_*.vtk"

        track_pattern =  "*-track_*.vtk"
        list = sorted(glob.glob(track_pattern))
        print "Looking with pattern " + track_pattern
        print "Found %d track files." % len(list)
        count = 0;
        for trackFile in list:

            # add plot
            # trackFile = conf['tracks_dir'] + os.path.sep + fname

            # plot the legend for the first one only
            if (count == 1) and legendFlag:
                trackPlotConf['PseudocolorAttributes']['legendFlag'] = 0
                # pp.pprint(trackPlotConf)

            # Plot the actual track data
            file = conf['tracks_dir'] + os.path.sep + trackFile
            print "Adding plot for " + file
            utils.addPseudolorPlot(file,trackPlotConf)

            count = count + 1

            # in case the script is being debugged, exit the script
            # after 10 tracks. This could be configured
            # if getValueForKeyPath(conf,'postprocessing.debugVisitScript') and count > 10

        # Restore flag value
        trackPlotConf['PseudocolorAttributes']['legendFlag'] = legendFlag
        # pp.pprint(trackPlotConf)

    print "Drawing plots"
    visit.DrawPlots()

    print "Saving image to %s" % os.getcwd()
    utils.saveImage("tracks",0)

    os.chdir(currentDirectory)

    return
