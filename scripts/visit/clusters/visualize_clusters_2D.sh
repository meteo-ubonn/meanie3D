#!/bin/bash

# ---------------------------------------------------------------
# This shell script constructs a python file from a template and
# executes the file in visit using the -s parameter.
#
# @author Jürgen Simon (juergen_simon@mac.com)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Prints usage and exits
# ---------------------------------------------------------------
usage() {
    echo "$0 <configuration file> <netcdf directory> <cluster directory> [resume]"
    echo "Creates a python script for visualizing 2D source and clusters for the meanie3D project."
    echo "    <configuration file> : json configuration file (check manual) (mandatory)"
    echo "    <netcdf directory> : directory containing the source data (mandatory)"
    echo "    <cluster directory> : directory containing the clustering results (mandatory)"
    echo "    [resume] : add this parameter to indicate that you wish to pick up where a previous"
    echo "               run left off. If this is not present, all previous data is erased"
    exit 1
}

# ---------------------------------------------------------------
# Checks for required environment variables and exits if
# they are not defined. The variables are:
#
# VISIT_EXECUTABLE : points to the 'visit' executable
# MEANIE3D_HOME : points to the checked out meanie3D directory
#                 which contains the 'visit' subdirectory
# ---------------------------------------------------------------
check_environ() {
    if [ "X${VISIT_EXECUTABLE}" = "X" ]
    then
        echo "Please set environment variable VISIT_EXECUTABLE"
        exit 1
    fi
    if [ "X${MEANIE3D_HOME}" = "X" ]
    then
        echo "Please set environment variable MEANIE3D_HOME"
        exit 1
    fi
}

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

# Check environment variables

check_environ

# mandatory arguments

if [ $# -lt 3 ]; then
    usage
fi

# optional arguments

RESUME="NO"
if [ "X$4" = "Xresume" ]
then
    RESUME="YES"
fi

# escape slashes

ESCAPED_CONFIG_FILE=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_NETCDF_DIR=$(echo $2 | sed -e "s/\//\\\\\//g")
ESCAPED_CLUSTER_DIR=$(echo $3 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

# construct python file from template

SCRIPTFILE="/tmp/cluster_visualization-$RANDOM.py"
cat $MEANIE3D_HOME/scripts/visit/clusters/visualize_clusters_2D.py \
| sed -e "s/P_CONFIGURATION_FILE/$ESCAPED_CONFIG_FILE/g" \
| sed -e "s/P_NETCDF_DIR/$ESCAPED_NETCDF_DIR/g" \
| sed -e "s/P_CLUSTER_DIR/$ESCAPED_CLUSTER_DIR/g" \
| sed -e "s/P_M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" \
| sed -e "s/P_RESUME/$RESUME/g" \
> ${SCRIPTFILE}

# run constructed file in visit

if [ "${RUN_VISIT_HEADLESS}" = "true" ]
then
    ${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPTFILE}
else
    ${VISIT_EXECUTABLE} -s ${SCRIPTFILE}
fi