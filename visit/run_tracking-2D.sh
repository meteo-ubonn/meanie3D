#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "run_tracking.sh <path to directory containing composite files>"
    echo "Creates a python script for cluster creation and tracking and runs it in Visit"
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

DL_PATH=$MEANIE3D_HOME/Release
SCRIPTFILE="/tmp/tracking-$RANDOM.py"
ESCAPED_SOURCE_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")
ESCAPED_DL_PATH=$(echo $DL_PATH | sed -e "s/\//\\\\\//g")
DETECTION_PARAMS="-r 12,12,200 --drf-threshold 0.5 -s 64 -t 20 -m 10"

cat $MEANIE3D_HOME/visit/run_tracking-2D.py | sed -e "s/SOURCE_DIR/$ESCAPED_SOURCE_DIR/g" | sed -e "s/DL_PATH/$ESCAPED_DL_PATH/g" | sed -e "s/M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" > $SCRIPTFILE
${VISIT_EXECUTABLE} -cli -nowin -s $SCRIPTFILE
#${VISIT_EXECUTABLE} -s $SCRIPTFILE