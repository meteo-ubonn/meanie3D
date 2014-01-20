#!/bin/bash

# Check environment variables first

if [ "X${VISIT_EXECUTABLE}" = "X" ]
then
    echo "Please set environment variable VISIT_EXECUTABLE"
    exit 0
fi
if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

# Check command line params

if [ "X$1" = "X" ]
then
    echo "visualize_tracks.py <dir> [<basename>]"
    echo "Creates a python script for cluster visualization in Visit and runs it."
    echo "The script finds all tracks in <dir> (looking for -track_*.vtk) and plots them."
    echo "Optional give basename, which will narrow to <basename>-track_*.vtk"
    exit 0
fi

# Escape params / variables

if [ "X$2" = "X" ]
then
    ESCAPED_BASENAME=$(echo $2 | sed -e "s/\//\\\\\//g")
else
    ESCAPED_BASENAME=""
fi
ESCAPED_SOURCE_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

# Create script

SCRIPT=$(echo "/tmp/visit_${RANDOM}.py")
echo ${SCRIPT}
cat ${MEANIE3D_HOME}/visit/visualize_tracks.py | sed -e "s/P_SOURCE_DIR/$ESCAPED_SOURCE_DIR/g" | sed -e "s/P_BASENAME/$ESCAPED_BASENAME/g" | sed -e "s/P_M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" > ${SCRIPT}

# Run script

echo "Starting Visit with script ${SCRIPT}"

# Headless
#${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPT}

# GUI
${VISIT_EXECUTABLE} -s ${SCRIPT}