#!/usr/bin/python

# ----------------------------------------------------------------------------
# run_tracking.py
#
# Python script for running a whole set of netcdf files through the
# clustering and tracking process.
# @author Juergen Simon (juergen_simon@mac.com)
# ----------------------------------------------------------------------------

# Import modules
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
# Prints usage and exits
# ----------------------------------------------------------------------------
def usage():
    print "run_clustering.py --configuration=<json file> --source=<netcdf directory> [--scale=<scale>] [--resume]"
    print "runs a complete set of netcdf files through the clustering/tracking"
    print "-c : json config file specifying variables etc."
    print "-f : directory containing the files. It is assumed that"
    print "           the files are in the correct order when sorted alphabetically."
    print "-s  : (optional) comma separated list of scale parameters. Overrides any scale values in the configuration."
    print "--resume,-r : if this flag is present, the algorithm assumes to resume"
    print "              operations where it last left off. If not present, previous"
    print "              results will be erased before starting"
    print "--help, -h  : print this message and exit."
    print "--version   : prints the version information and exits"
    sys.exit(1)
    return

# ----------------------------------------------------------------------------
# Prints version info and exits
# ----------------------------------------------------------------------------

def print_version():
    DYLD_LIBRARY_PATH="/usr/local/lib:/usr/lib"
    bin_prefix    = "export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:"+DYLD_LIBRARY_PATH+";"
    print "run_clustering.py : "
    print meanie3D.get_version()
    print "meanie3D-detect   : "
    call( bin_prefix + "/usr/local/bin/" + "meanie3D-detect --version", shell=True)
    print "meanie3D-track    : "
    call( bin_prefix + "/usr/local/bin/" + "meanie3D-track --version", shell=True)
    sys.exit(1)
    return

# ----------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------
def main(argv):
    
    # Parse command line

    try:
        opts, args = getopt.getopt(argv, "c:f:s:o:rh", ["resume","help","version"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    scales = []
    resume = False
    config_file = ""
    output_dir = "."
    netcdf_dir = ""
    num_params = 0
    
    for o, a in opts:
        
        print o
        print a
        print
        
        if o == "-c":
            config_file = a
            num_params = num_params + 1
        
        elif o == "-f":
            netcdf_dir = a
            num_params = num_params + 1
        
        elif o == "-s":
            scales = str(a).split(',')

        elif o == "-o":
            output_dir = a
        
        elif o in ["--resume","-r"]:
            resume = True
        
        elif o in ["--help"]:
            usage()
            sys.exit()

        elif o in ["--version"]:
            print_version()

        else:
            usage()
                
    if num_params < 2:
        print num_params
        usage()

    # Parse configuration data and expand

    configuration = meanie3D.load_configuration(config_file);
    configuration["NETCDF_DIR"] = netcdf_dir
    configuration["OUTPUT_DIR"] = output_dir
    configuration["M3D_HOME"] = MEANIE3D_HOME
    configuration["RESUME"] = resume

    #
    # run the actual clustering/tracking script
    #

    if not scales:
        meanie3D.run_tracking(configuration,-1)

    else:
        for scale in scales:
            configuration["SCALE"] = scale
            meanie3D.run_tracking(configuration,-1)

    return

# ----------------------------------------------------------------------------
# entry point
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    main(sys.argv[1:])