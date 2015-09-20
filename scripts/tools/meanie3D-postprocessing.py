#!/usr/bin/python

## meanie3D-tracking.py
#
# Python script for running a whole set of netcdf files through the clustering and tracking process.
# \author Juergen Simon (juergen.simon@uni-bonn.de)

import glob
import os
import sys
import getopt
from subprocess import call

# dynamically append search path
# to include the meanie3D modules
# Parameters
MEANIE3D_HOME = os.getenv("MEANIE3D_HOME","NOT_SET");
if MEANIE3D_HOME == "NOT_SET":
    print "ERROR: environment variable MEANIE3D_HOME is not set."
    sys.exit(2)
sys.path.append(MEANIE3D_HOME+"/scripts/python-modules")
import meanie3D

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
    print "meanie3D.py:"
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
def run_trackstats(configuration,directory):
    print "Running meanie3D-trackstats for %s" % directory
    return

# ----------------------------------------------------------------------------
## Runs gnuplot to produce .eps files
# \param configuration
# \param directory
def plot_trackstats(configuration,directory):
    print "Plotting .eps files for %s" % directory
    return

# ----------------------------------------------------------------------------
## Creates a python script for Visit and runs Visit with the script.
# Visualises the tracks found.
# \param configuration
# \param directory
def visualise_tracks(configuration,directory):
    print "Visualising tracks for %s" % directory
    return

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
    configuration["netcdf_dir"] = netcdf_dir
    configuration["m3d_home"] = MEANIE3D_HOME

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
            run_trackstats(configuration, directory)

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