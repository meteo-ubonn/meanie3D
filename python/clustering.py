#!/usr/bin/python

import argparse
import os
import sys

# Check environment variables

MEANIE3D_HOME=os.getenv("MEANIE3D_HOME")
if MEANIE3D_HOME=="none":
    print "Please set MEANIE3D_HOME"
    exit(-1)

VISIT_EXECUTABLE=os.getenv("VISIT_EXECUTABLE")
if VISIT_EXECUTABLE=="none":
    print "Please set VISIT_EXECUTABLE"
    exit(-1)

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def main(argv):

    # processing 
    source_dir="."
    dest_dir="."
    clustering_options=""
    dimensions=""
    vtk_dimensions=""
    variables=""
    lower_thresholds=""
    upper_thresholds=""
    

    # visualization
    visualize_source=False
    visualize_scalespace=False
    visualize_weight_function=False
    visualize_meanshift=False
    ancillary_data_file=""
    ancillary_data_vars=[];

    ap = argparse.ArgumentParser()
    ap.add_argument("--source",type=string,default=".",help="Directory containing the NetCDF files to cluster")
    ap.add_argument("--destination",type=string,default=".",help="Directory to write the resulting cluster files to")
    
    # iterate over the raw argument list and stick
    # everything that is not 


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


if __name__ == "__main__":
    main(sys.argv[1:])
