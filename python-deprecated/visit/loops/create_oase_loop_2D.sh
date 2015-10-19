#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "visualize_clusters.sh <netcdf directory> <variable_name>"
    echo "Creates a python script for visualizing variables in 2D OASE composite files"
    exit 0
fi

if [ "X$2" = "X" ]
then
    echo "visualize_clusters.sh <netcdf directory> <variable_name>"
    echo "Creates a python script for visualizing variables in 2D OASE composite files"
    exit 0
fi

COLOR_TABLE="hot_desaturated"
INVERT_COLOR_TABLE="0"

echo "NUMBER OF PARAMS $#"

if [ $# -gt 2 ]
then
    echo "COLOR_TABLE is set to " + $3
    COLOR_TABLE=$3
else
    echo "COLOR_TABLE is set to default " + $COLOR_TABLE
fi

if [ $# -gt 3 ]
then
    echo "INVERT_COLOR_TABLE is set to " + $4
    INVERT_COLOR_TABLE=$4
else
    echo "INVERT_COLOR_TABLE is set to default " + $INVERT_COLOR_TABLE
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

SCRIPTFILE="/tmp/cluster_visualization-$RANDOM.py"
ESCAPED_NETCDF_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/scripts/visit/loops/create_oase_loop_2D.py | sed -e "s/SOURCE_DIR_P/$ESCAPED_NETCDF_DIR/g" | sed -e "s/VARIABLE_P/$2/g" | sed -e "s/COLOR_TABLE_P/$COLOR_TABLE/g" | sed -e "s/COLOR_TABLE_INVERT_P/$INVERT_COLOR_TABLE/g" | sed -e "s/MEANIE3D_HOME_P/$ESCAPED_MEANIE3D_HOME/g" > ${SCRIPTFILE}

if [ "${RUN_VISIT_HEADLESS}" = "true" ]
then
    ${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPTFILE}
else
    ${VISIT_EXECUTABLE} -s ${SCRIPTFILE}
fi
