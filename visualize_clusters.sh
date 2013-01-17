#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "visualize_clusters.sh <path to source file>"
    echo "Creates a python script for cluster visualization in Visit and runs it"
    exit 0
fi

if [ "X${VISIT_EXECUTABLE}" = "X" ]
then
    echo "Please set environment variable VISIT_EXECUTABLE"
    exit 0
fi

ESCAPED_SOURCE_FILE=$(echo $1 | sed -e "s/\//\\\\\//g")
cat cluster_visualization.template | sed -e "s/SOURCE_FILE/$ESCAPED_SOURCE_FILE/g" > visit_script.py
${VISIT_EXECUTABLE} -s visit_script.py