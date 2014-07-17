#!/bin/bash

if [ "X${VISIT_EXECUTABLE}" = "X" ]; then
    echo "Please set environment variable VISIT_EXECUTABLE"
    exit 0
fi

if [ "X${MEANIE3D_HOME}" = "X" ]; then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

OVERLAY="NO"
RESUME_AT_FILE=""

if [ "X$1" = "X" ]; then
    echo "visualize_clusters.sh <netcdf directory> <cluster directory> [<resume at file>] [overlay]"
    echo "Creates a python script for cluster visualization in Visit and runs it"
    exit 0
fi

if [ "X$2" = "X" ]; then
    echo "visualize_clusters.sh <netcdf directory> <cluster directory> [<resume at file>] [overlay]"
    echo "Creates a python script for cluster visualization in Visit and runs it"
    exit 0
fi

if [ "$3" = "overlay" ]; then
    OVERLAY="YES"
else
    RESUME_AT_FILE="$3"
fi

if [ "$4" = "overlay" ]; then
    OVERLAY="YES"
fi

echo "Resuming at file ${RESUME_AT_FILE}"
echo "Producing overlays (shell): ${OVERLAY}"

SCRIPTFILE="/tmp/cluster_visualization-$RANDOM.py"
ESCAPED_NETCDF_DIR=$(echo $1 | sed -e "s/\//\\\\\//g")
ESCAPED_CLUSTER_DIR=$(echo $2 | sed -e "s/\//\\\\\//g")
ESCAPED_MEANIE3D_HOME=$(echo $MEANIE3D_HOME | sed -e "s/\//\\\\\//g")
ESCAPED_RESUME_FILE=$(echo $RESUME_AT_FILE | sed -e "s/\//\\\\\//g")

cat $MEANIE3D_HOME/visit/visualize_results/visualize_oase_clusters_2D.py | sed -e "s/P_NETCDF_DIR/$ESCAPED_NETCDF_DIR/g" | sed -e "s/P_CLUSTER_DIR/$ESCAPED_CLUSTER_DIR/g" | sed -e "s/P_M3D_HOME/$ESCAPED_MEANIE3D_HOME/g" | sed -e "s/P_OVERLAY/$OVERLAY/g" | sed -e "s/P_RESUME_AT_FILE/${ESCAPED_RESUME_FILE}/g" > ${SCRIPTFILE}

if [ "${RUN_VISIT_HEADLESS}" = "true" ]
then
    ${VISIT_EXECUTABLE} -cli -nowin -s ${SCRIPTFILE}
else
    ${VISIT_EXECUTABLE} -s ${SCRIPTFILE}
fi

