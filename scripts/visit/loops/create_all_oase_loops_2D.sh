#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "create_all_oase_loops_2D.sh <netcdf directory>"
    echo "Visualize all variables in 2D OASE composite files"
    exit 0
fi

if [ "X${VISIT_EXECUTABLE}" = "X" ]
then
    echo "Please set environment variable VISIT_EXECUTABLE"/Users/simon/Desktop/tracking_0144-cropped.png
    exit 0
fi

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

SCRIPTFILE="/tmp/cluster_visualization-$RANDOM.py"
ESCAPED_NETCDF_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/scripts/visit/loops/create_all_oase_loops_2D.py | sed -e "s/SOURCE_DIR_P/$ESCAPED_NETCDF_DIR/g" | sed -e "s/MEANIE3D_HOME_P/$ESCAPED_MEANIE3D_HOME/g" > ${SCRIPTFILE}

if [ "${RUN_VISIT_HEADLESS}" = "true" ]
then
    ${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPTFILE}
else
    ${VISIT_EXECUTABLE} -s ${SCRIPTFILE}
fi
