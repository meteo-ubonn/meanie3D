#!/usr/bin/python

## meanie3D-tracking.py
#
# Python script for running a whole set of netcdf files through the clustering and tracking process.
# \author Juergen Simon (juergen.simon@uni-bonn.de)

# System packages
import glob
import os
import sys
import visit

# Own packages
from . import __have_visit__
from . import __visitPath__

from meanie3D.app import utils
from meanie3D.app import external

# make sure external commands are available
external.find_ext_cmds(['meanie3D-cfm2vtk','meanie3D-trackstats','gnuplot'])

# ----------------------------------------------------------------------------
## Checks the consistency of the given configuration and hotfixes
# problems it finds.
# \param configuration
#
def check_configuration(configuration):

    # Check the configuration
    dimensions = utils.getValueForKeyPath(configuration,'data.dimensions')
    if (not dimensions):
        print "ERROR:configuration must contain 'dimensions'"
        return -1
    if (not len(dimensions) in [2,3]):
        print "ERROR:Can only process 2D or 3D"
        return -1

    # Check that visualiseTracks has .vtk to work from
    if utils.getValueForKeyPath(configuration,'postprocessing.tracks.visualiseTracks') and not utils.getValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.vtk_tracks'):
        print "WARNING: tracks.visualiseTracks = True but tracks.meanie3D-trackstats.vtk_tracks = False. Correcting."
        utils.setValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.vtk_tracks',True)

    # Complement vtk_dimensions to make our life a little easier down the road.
    vtkDimensions = utils.getValueForKeyPath(configuration,'data.vtkDimensions')
    if vtkDimensions:
        vtkDimString = ",".join(vtkDimensions)
        utils.setValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.vtkDimensions',vtkDimString)

    # Make sure that plotStats has gnuplot files to work with
    if utils.getValueForKeyPath(configuration,'postprocessing.tracks.plotStats') and not utils.getValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.gnuplot'):
        print "WARNING: track.plot_stats = True but tracks.meanie3D-trackstats.gnuplot = False. Correcting."
        utils.setValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.gnuplot',True)

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
## Runs the meanie3D-trackstats command.
# \param configuration
# \param directory
# \return True if the stats were created, False else
def run_trackstats(configuration,directory):
    pconf = configuration['postprocessing']
    if not 'tracks' in pconf:
        return False
    tracks = pconf['tracks']
    if not 'meanie3D-trackstats' in tracks:
        return False

    print "Running meanie3D-trackstats for %s" % directory
    conf = tracks['meanie3D-trackstats']

    # Assemble command line params
    params = []
    if (conf['dictionary'] == True):
        params.append("-t")
    if (conf['gnuplot'] == True):
        params.append("-g")
    if (conf['length'] == True):
        params.append("--create-length-statistics")
    if (conf['length_classes']):
        params.append("--length-histogram-classes=%s" % conf['length_classes'])
    if (conf['speed'] == True):
        params.append("--create-speed-statistics")
    if (conf['speed_classes']):
        params.append("--speed-histogram-classes=%s" % conf['speed_classes'])
    if (conf['direction'] == True):
        params.append("--create-direction-statistics")
    if (conf['direction_classes']):
        params.append("--direction-histogram-classes=%s" % conf['direction_classes'])
    if (conf['size'] == True):
        params.append("--create-cluster-statistics")
    if (conf['size_classes']):
        params.append("--cluster-histogram-classes=%s" % conf['size_classes'])
    if (conf['cumulated'] == True):
        params.append("--create-cumulated-size-statistics")
    if (conf['cumulated_classes']):
        params.append("--size-histogram-classes=%s" % conf['cumulated_classes'])
    if (conf['vtk_tracks'] == True):
        params.append("--write-center-tracks-as-vtk ")
    if (conf['vtkDimensions']):
        params.append("--vtk-dimensions=%s" % conf['vtkDimensions'])

    params.append("-s netcdf")

    return_code = -1
    os.chdir(directory)
    try:
        return_code = external.execute_command("meanie3D-trackstats"," ".join(params))
    except:
        print "ERROR:%s" % sys.exc_info()[0]
        raise

    os.chdir("..")

    return (return_code == 0)

# ----------------------------------------------------------------------------
## Runs gnuplot to produce .eps files
# \param configuration
# \param directory
def plot_trackstats(configuration,directory):
    pconf = configuration['postprocessing']
    if not 'tracks' in pconf:
        return False
    tracks = pconf['tracks']
    if not 'meanie3D-trackstats' in tracks:
        return False

    conf = tracks['meanie3D-trackstats']
    print "Plotting .eps files for %s" % directory

    os.chdir(directory)

    # create a tempfile for gnuplot
    f = open('plot_stats.gp','w')

    # track length
    if (conf['length'] == True):
        f.write('set term postscript\n')
        f.write('set output "lengths-hist.eps"\n')
        f.write('set title "Distribution of track length"\n')
        f.write('set xlabel "Track length in [#steps]"\n')
        f.write('set ylabel "Number of tracks"\n')
        f.write('set xtics 5\n')
        f.write('plot "lengths-hist.txt" with boxes\n')

    # sizes
    if (conf['size']):
        f.write('# cluster size\n')
        f.write('set output "sizes-hist.eps"\n')
        f.write('set title "Distribution of cluster size"\n')
        f.write('set xlabel "log(cluster size in [#gridpoints])"\n')
        f.write('set ylabel "log(number of clusters)"\n')
        f.write('set logscale y\n')
        f.write('set logscale x\n')
        f.write('set xtics auto\n')
        f.write('plot "sizes-hist.txt" with boxes\n')
        f.write('unset logscale y\n')
        f.write('unset logscale x\n')

    # speed
    if (conf['speed']):
        f.write('set output "speeds-hist.eps"\n')
        f.write('set title "Distribution of cluster speeds"\n')
        f.write('set xlabel "Cluster speed in [m/s]"\n')
        f.write('set ylabel "Number of clusters"\n')
        f.write('set xtics 2\n')
        f.write('plot "speeds-hist.txt" with boxes\n')

    # directions
    if (conf['direction']):
        f.write('set output "directions-hist.eps"\n')
        f.write('set title "Distribution of tracking direction"\n')
        f.write('set xlabel "Cluster direction in [deg]"\n')
        f.write('set ylabel "Number of clusters"\n')
        f.write('set xtics 15\n')
        f.write('plot "directions-hist.txt" with boxes\n')

    f.close()

    return_code = -1
    try:
        return_code = external.execute_command("gnuplot","plot_stats.gp")
    except:
        print "ERROR:%s" % sys.exc_info()[0]

    os.chdir("..")
    return (return_code == 0)

# ----------------------------------------------------------------------------
## Creates a python script for Visit and runs Visit with the script.
# Visualises the tracks found.
# \param configuration
# \param directory
def visualise_tracks(configuration,directory):
    print "Visualising tracks for %s" % directory
    if not __have_visit__:
        return

    try:
        visit.Launch(vdir=__visitPath__)
    except:

    # # Find template
    # templatePath = os.path.join(os.path.split(__file__)[0], "templates/tracks_visit.py")
    # replacements = {
    #     'P_TRACKS_DIR' : os.path.abspath(directory),
    #     'P_M3D_HOME' : meanie3D.__file__,
    #     'P_CONFIGURATION_FILE' : os.path.abspath(configuration['config_file'])
    # }
    #
    # scriptFilename = tempfile.mktemp() + ".py"
    # print "\tWriting python script for visualisation to: " + scriptFilename
    # with open(templatePath) as infile, open(scriptFilename, 'w') as outfile:
    #     for line in infile:
    #         for src, target in replacements.iteritems():
    #             line = line.replace(src, target)
    #         outfile.write(line)
    # print "\tDone."
    #
    # # Compile command line params for visualisation
    # params = "-s %s" % scriptFilename
    # runHeadless = utils.getValueForKeyPath(configuration,'postprocessing.runVisitHeadless')
    # if runHeadless:
    #     params = "-cli -nowin %s" % params
    #
    # # change into directory
    # os.chdir(directory)
    # # Run visualisation
    # returnCode = -1
    # try:
    #     returnCode = external.execute_command("visualisation",params)
    # except:
    #     print "ERROR:%s" % sys.exc_info()[0]
    #
    # # Change back out
    # os.chdir('..')
    # return (returnCode == 0)

# ----------------------------------------------------------------------------
## Runs cross-scale analysis
# \param configuration
# \param directory
def run_comparison(configuration):
    print "Running cross-scale comparison"
    return

# ----------------------------------------------------------------------------
## Creates a python script for Visit and runs Visit with the script.
# Visualises the clusters, creates movies etc.
# \param configuration
# \param directory
def visualise_clusters(configuration,directory):
    print "Visualising clusters for %s" % directory
    return

# ----------------------------------------------------------------------------
## Removes any files that are not part of the results and were produced
# in any of the subsequent steps. C
# \param configuration
# \param directory
def remove_junk(configuration,directory):
    print "Cleaning temporary files for %s" % directory
    os.chdir(directory)

    if utils.getValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.vtk_tracks'):
        for filename in glob.glob('*.vtk') :
            os.remove(filename)

    if utils.getValueForKeyPath(configuration,'postprocessing.tracks.meanie3D-trackstats.gnuplot'):
        for filename in glob.glob('*.gp') :
            os.remove(filename)

    os.chdir("..")
    return

##
# Runs postprocessing steps according to the section 'postprocessing' in the
# given configuration.
# \param:directory
# \param:configuration
#
def run(configuration):

    # Check configuration section 'postprocessing'
    check_configuration(configuration)

    if not 'postprocessing' in configuration:
        return



    # In case scale parameters were given, the output dirs are scaleXYZ.
    # Otherwise it's 'clustering'. To be safe, iterate over both
    directories = sorted(glob.glob("scale*"),reverse=True)
    if (os.path.isdir("clustering")):
        directories.append("clustering")

    print "Processing directories: %s" % str(directories)
    for directory in directories:
        if os.path.isdir(directory):
            print "Processing %s" % directory
        else:
            print "ERROR:%s is not a directory!" % directory
            continue

        # run the track statistics
        if (run_trackstats(configuration, directory)):

            # run the stats plotting
            if utils.getValueForKeyPath(configuration,'postprocessing.tracks.plotStats'):
                plot_trackstats(configuration,directory);

            # run the track visualisations
            if utils.getValueForKeyPath(configuration, 'postprocessing.tracks.visualiseTracks'):
                visualise_tracks(configuration, directory)

        if utils.getValueForKeyPath(configuration,'postprocessing.clusters.visualiseClusters'):
            visualise_clusters(configuration,directory)

        remove_junk(configuration,directory)

    if (utils.getValueForKeyPath(configuration,'postprocessing.tracks.runScaleComparison')):
        run_comparison(configuration)

# ----------------------------------------------------------------------------

