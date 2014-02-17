#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "run_tracking.sh <path to directory containing composite files> <range[m]>"
    echo "Creates a python script for cluster creation and tracking and runs it in Visit"
    exit 0
fi

if [ "X$2" = "X" ]
then
    echo "run_tracking.sh <path to directory containing composite files> <range[m]>"
    echo "Creates a python script for cluster creation and tracking and runs it in Visit"
    exit 0
exit 0

fi

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

SCRIPTFILE="/tmp/tracking-$RANDOM.py"
ESCAPED_SOURCE_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/visit/tracking/run_rico_tracking_2D.py | sed -e "s/PARAM_SOURCE_DIR/$ESCAPED_SOURCE_DIR/g" | sed -e "s/M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" | sed -e "s/PARAM_SCALE/$2/g" > $SCRIPTFILE

python ${SCRIPTFILE} 