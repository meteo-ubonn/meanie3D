#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "create_radolan_loop.sh <path to directory containing composite files>"
    echo "Creates a python script for visualizing all scans in the given directory and runs it in visit"
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

#DL_PATH=$MEANIE3D_HOME/Release
DL_PATH=/usr/local/lib

SCRIPTFILE="/tmp/tracking-$RANDOM.py"
ESCAPED_SOURCE_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")
ESCAPED_DL_PATH=$(echo $DL_PATH | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/visit/create_radolan_loop.py | sed -e "s/SOURCE_DIR/$ESCAPED_SOURCE_DIR/g" | sed -e "s/DL_PATH/$ESCAPED_DL_PATH/g" | sed -e "s/M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" > $SCRIPTFILE
${VISIT_EXECUTABLE} -cli -nowin -s $SCRIPTFILE
