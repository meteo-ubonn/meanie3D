#!/usr/bin/python

## meanie3D-tracking.py
#
# Python script for running a whole set of netcdf files through the clustering and tracking process.
# \author Juergen Simon (juergen.simon@uni-bonn.de)

import glob
import os
import sys
import getopt
import tempfile
from subprocess import call

# dynamically append search path
# to include the meanie3D modules
# Parameters
MEANIE3D_HOME = os.getenv("MEANIE3D_HOME","NOT_SET");
if MEANIE3D_HOME == "NOT_SET":
    print "ERROR: environment variable MEANIE3D_HOME is not set."
    sys.exit(2)
sys.path.append(MEANIE3D_HOME+"/python/python-modules")
from meanie3D import external
import meanie3D

# make sure external commands are available
external.find_ext_cmds(['meanie3D-cfm2vtk','meanie3D-trackstats','gnuplot','visit'])

# ----------------------------------------------------------------------------
## Prints usage and exits
#
def usage():
    print "meanie3D-postprocessing.py -c <json configuration file> -s <netcdf directory> [--help,-h] [--version]"
    print "Runs a number of postprocessing steps on the results of meanie3D-tracking.py"
    print "-c : json config file specifying variables etc."
    print "-s : directory containing the NetCDF source data (see meanie3D-tracking.py --help)"
    print "--json-example : prints an example .json configuration and exits"
    print "--help, -h  : print this message and exit."
    print "--version   : prints the version information and exits"
    sys.exit(1)
    return
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
## Prints the configuration file .json file description and exits
#
def print_configuration_format():
    print '''{'''
    print '''    "description" : "Description of file content",'''
    print '''    "tracks" : {'''
    print '''       "dictionary" : true,        /* See --write-track-dictionary on meanie3D-trackstats */'''
    print '''       "gnuplot" : true,           /* true|false. See --write-gnuplot-files on meanie3D-trackstats */'''
    print '''       "length" : true,            /* true|false. See --create-length-statistics on meanie3D-trackstats */'''
    print '''       "length_classes" : null,    /* null or []. See --length-histogram-classes on meanie3D-trackstats */'''
    print '''       "speed" : true,             /* --create-speed-statistics  on meanie3D-trackstats */'''
    print '''       "speed_classes" : null,     /* null or []. See -speed-histogram-classes  on meanie3D-trackstats */'''
    print '''       "direction": true,          /* true|false. See --create-direction-statistics on meanie3D-trackstats */'''
    print '''       "direction_classes" : null, /* null or []. See --direction-histogram-classes  on meanie3D-trackstats */'''
    print '''       "size": true,               /* true|false. See --create-cluster-statistics on meanie3D-trackstats */'''
    print '''       "size_classes" : null,      /* null or []. See --cluster-histogram-classes on meanie3D-trackstats */'''
    print '''       "cumulated" : false,        /* true|false See --create-cumulated-size-statistics on meanie3D-trackstats */'''
    print '''       "cumulated_classes" : null, /* null or []. See --size-histogram-classes arg on meanie3D-trackstats */'''
    print '''       "vtk_tracks" : true,        /* true|false. See --write-center-tracks-as-vtk on meanie3D-trackstats */'''
    print '''       "vtk_dimensions" : "x,y",   /* null or "x,y,..". See --vtk-dimensions on meanie3D-trackstats */'''
    print '''       "plot_stats" : true,        /* true|false. Plots .eps files of the stats. Implies "gnuplot" : true*/'''
    print '''       "visualise_tracks" : true,  /* If true, all tracks are plotted on one image. Implies"vtk_tracks": true.*/'''
    print '''       "scale_comparison" : true   /* If true, available stats are compared across scales. */'''
    print '''    },'''
    print '''    "clusters" :'''
    print '''    {'''
    print '''       "variables":["var1","var2"],         /* Comma-separated list of variables to plot on the source movies. They are plotted in that order.*/'''
    print '''       "lower_tresholds" : [20,50],         /* Comma-separated list of lower thresholds for each variable.*/'''
    print '''       "upper_tresholds" : [120.0,100],     /* Comma-separated list of upper thresholds for each variable.*/'''
    print '''       "var_min" : [0,0],                   /* Comma-separated list of minimum values for the color scale for each variable.*/'''
    print '''       "var_max" : [75.0,100.0],            /* Comma-separated list of maximum values for the color scale for each variable.*/'''
    print '''       "with_background_gradient" : true,   /* true|false. If true, a background gradient is plotted on each image.*/'''
    print '''       "with_topography" : true,            /* true|false. If true, the country's topography from the map data file is plotted as background.*/'''
    print '''       "with_rivers_and_boundaries" : true, /* true|false. If true, the country's boundaries and rivers are plotted as background.*/'''
    print '''       "with_source_backround" : true,      /* true|false. If true, the source data is plotted before the clusters are added.*/'''
    print '''       "with_datetime" :  true,             /* true|false. If true. a timestamp is added in the top right corner. */'''
    print '''       "create_source_movie" : true,        /* true|false. If true, a movie is created from the source data.*/'''
    print '''       "create_clusters_movie" : true,      /* true|false. If true, a movie is created from the cluster data.*/'''
    print '''       "movie_formats": ["gif","m4v"],      /* Comma-separated lists of file extensions for the movies. Extension must be understood by imagemagick's convert command. */'''
    print '''       "movie_convert_params" : ["",""],    /* parameters to imagemagick's convert commamnd to use for each movie format parameter. */'''
    print '''       "grid_extent" : "national",          /* Denotes the type of grid to use (radolan is "national", Oase's Bonn-Juelich format is "local"). */'''
    print '''       "meanie3D-cfm2vtk" : "-t cluster -v RX --vtk-dimensions x,y",   /* Parameters to pass to meanie3D-cfm2vtk.*/'''
    print '''       "colortables" : ["hot","blues"],     /* Visit colortable names to use for each variable. */'''
    print '''       "colortables_invert_flags" : [0,1], /* Array of flags indicating if colortable should be inverted for each variable.*/'''
    print '''       "opacity" : [1.0,0.33]              /* Opacity for each variable. Example: [1.0]. */'''
    print '''    }'''
    print '''}'''
    sys.exit(1)
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
## Prints version info and exits
#
def print_version():
    DYLD_LIBRARY_PATH="/usr/local/lib:/usr/lib"
    bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
    print "utils.py:"
    print meanie3D.get_version()
    print "meanie3D-trackstats:"
    call( bin_prefix + "/usr/local/bin/" + "meanie3D-trackstats --version", shell=True)
    print "meanie3D-cfm2vtk:"
    call( bin_prefix + "/usr/local/bin/" + "meanie3D-cfm2vtk --version", shell=True)
    sys.exit(1)
    return
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
## Checks the consistency of the given configuration and hotfixes
# problems it finds.
# \param configuration
#
def check_configuration(configuration):

    if (configuration['tracks']):

        # visualise_tracks -> vtk_tracks:true
        if (configuration['tracks']["visualise_tracks"] and not configuration['tracks']['vtk_tracks']):
            print "WARNING: tracks.visualise_tracks = True -> tracks.vtk_tracks = True"
            configuration['tracks']['vtk_tracks'] = True

        # visualise_tracks -> vtk_tracks:true
        if (configuration['tracks']["plot_stats"] and not configuration['tracks']['gnuplot']):
            print "WARNING: track.plot_stats = True -> track.gnuplot = True"
            configuration['tracks']['gnuplot'] = True

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
## Runs the meanie3D-trackstats command.
# \param configuration
# \param directory
# \return True if the stats were created, False else
def run_trackstats(configuration,directory):
    print "Running meanie3D-trackstats for %s" % directory
    conf = configuration['tracks']
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
    if (conf['vtk_dimensions']):
        params.append("--vtk-dimensions=%s" % conf['vtk_dimensions'])

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
    print "Plotting .eps files for %s" % directory
    os.chdir(directory)

    # create a tempfile for gnuplot
    f = open('plot_stats.gp','w')

    # track length
    if (configuration['tracks']['length'] == True):
        f.write('set term postscript\n')
        f.write('set output "lengths-hist.eps"\n')
        f.write('set title "Distribution of track length"\n')
        f.write('set xlabel "Track length in [#steps]"\n')
        f.write('set ylabel "Number of tracks"\n')
        f.write('set xtics 5\n')
        f.write('plot "lengths-hist.txt" with boxes\n')

    # sizes
    if (configuration['tracks']['size']):
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
    if (configuration['tracks']['speed']):
        f.write('set output "speeds-hist.eps"\n')
        f.write('set title "Distribution of cluster speeds"\n')
        f.write('set xlabel "Cluster speed in [m/s]"\n')
        f.write('set ylabel "Number of clusters"\n')
        f.write('set xtics 2\n')
        f.write('plot "speeds-hist.txt" with boxes\n')

    # directions
    if (configuration['tracks']['direction']):
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
        raise

    os.chdir("..")
    return (return_code == 0)

# ----------------------------------------------------------------------------
## Creates a python script for Visit and runs Visit with the script.
# Visualises the tracks found.
# \param configuration
# \param directory
def visualise_tracks(configuration,directory):
    print "Visualising tracks for %s" % directory

    # Find template
    templatePath = os.path.expandvars("${MEANIE3D_HOME}/python/tools/templates/visualise_tracks.py")
    if not os.path.exists(templatePath):
        print "ERROR: could not find script template 'visualise_tracks.py' in $MEANIE3D_HOME"
        return -1

    replacements = {
        'P_TRACKS_DIR' : os.path.abspath(directory),
        'P_M3D_HOME' : configuration['m3d_home'],
        'P_CONFIGURATION_FILE' : os.path.abspath(configuration['config_file'])
    }

    scriptFilename = tempfile.mktemp() + ".py"
    print "\tWriting python script for visit to: " + scriptFilename
    with open(templatePath) as infile, open(scriptFilename, 'w') as outfile:
        for line in infile:
            for src, target in replacements.iteritems():
                line = line.replace(src, target)
            outfile.write(line)
    print "\tDone."

    # Compile command line params for visit
    params = "-s %s" % scriptFilename
    debugVisitScript = configuration['tracks']['debugVisitScript']
    runHeadless = bool(os.path.expandvars("$RUN_VISIT_HEADLESS")) and not debugVisitScript
    if runHeadless:
        params = "-cli -nowin %s" % params

    # change into directory
    os.chdir(directory)
    # Run visit
    external.execute_command("visit",params)
    # Change back out
    os.chdir('..')

    return 0

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

    if (configuration['tracks']['vtk_tracks']):
        for filename in glob.glob('*.vtk') :
            os.remove(filename)

    if (configuration['tracks']['gnuplot']):
        for filename in glob.glob('*.gp') :
            os.remove(filename)

    os.chdir("..")
    return

# ----------------------------------------------------------------------------
## Main function
def main(argv):
    # Parse command line
    try:
        opts, args = getopt.getopt(argv, "c:s:h", ["json-example","help","version"])
    except getopt.GetoptError as detail:
        print detail
        sys.exit(2)

    config_file = ""
    netcdf_dir = ""
    num_params = 0

    for o, a in opts:
        if o == "-c":
            config_file = a
            num_params = num_params + 1
        
        elif o == "-s":
            netcdf_dir = a
            num_params = num_params + 1

        elif o in ["--json-example"]:
            print_configuration_format()

        elif o in ["--help"]:
            usage()

        elif o in ["--version"]:
            print_version()

        else:
            usage()
                
    if num_params < 2:
        print num_params
        usage()

    # Parse configuration data and expand
    configuration = meanie3D.load_configuration(config_file);
    configuration['netcdf_dir'] = netcdf_dir
    configuration['m3d_home'] = os.path.expandvars("${MEANIE3D_HOME}")
    configuration['config_file'] = config_file

    # Check configuration consistency
    check_configuration(configuration)

    # In case scale parameters were given, the output dirs are scaleXYZ.
    # Otherwise it's 'clustering'. To be safe, iterate over both

    directories = sorted(glob.glob("scale*"))
    if (os.path.isdir("clustering")):
        directories.append("clustering")

    print "Processing directories: %s" % str(directories)
    for directory in directories:
        if os.path.isdir(directory):
            print "Processing %s" % directory
        else:
            print "ERROR:%s is not a directory!" % directory
            continue

        if (configuration['tracks']):

            # run the track statistics
            if (run_trackstats(configuration, directory)):

                # run the stats plotting
                if (configuration['tracks']['gnuplot']):
                    plot_trackstats(configuration,directory);

                # run the track visualisations
                if (configuration['tracks']['visualise_tracks']):
                    visualise_tracks(configuration, directory)

        if (configuration['clusters']):
            visualise_clusters(configuration,directory)

        remove_junk(configuration,directory)

    if (configuration['tracks'] and configuration['tracks']['scale_comparison']):
        run_comparison(configuration)

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# entry point
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    main(sys.argv[1:])