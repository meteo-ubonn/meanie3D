#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "visualize_tb4_clusters.sh <source file> <cluster file> <variable>"
    echo "Creates a python script for cluster visualization in Visit and runs it"
    exit 0
fi

if [ "X$2" = "X" ]
then
    echo "visualize_tb4_clusters.sh <source file> <cluster file> <variable>"
    echo "Creates a python script for cluster visualization in Visit and runs it"
    exit 0
fi

if [ "X$3" = "X" ]
then
    echo "visualize_tb4_clusters.sh <source file> <cluster file> <variable>"
    echo "Creates a python script for cluster visualization in Visit and runs it"
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

SCRIPTFILE="/tmp/tb4_visualization-$RANDOM.py"
ESCAPED_NETCDF_FILE=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_CLUSTER_FILE=$(echo $2 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/visit/tb4/visualize_tb4_2D.py | sed -e "s/P_NETCDF_FILE/$ESCAPED_NETCDF_FILE/g" | sed -e "s/P_CLUSTER_FILE/$ESCAPED_CLUSTER_FILE/g" | sed -e "s/P_VAR_NAME/$3/g" | sed -e "s/P_M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" > ${SCRIPTFILE}

if [ "${RUN_VISIT_HEADLESS}" = "true" ]
then
    ${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPTFILE}
else
    ${VISIT_EXECUTABLE} -s ${SCRIPTFILE}
fi