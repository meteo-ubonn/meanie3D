#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "trackstats-conrad.sh <conrad short format cell file>"
    exit 0
fi


if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
exit 0
fi

meanie3D-trackstats-conrad --write-track-dictionary --write-gnuplot-files \
    --create-length-statistics --create-speed-statistics \
    --create-direction-statistics --create-cluster-statistics \
    --write-center-tracks-as-vtk \
    -f $1

# plot stats
${MEANIE3D_HOME}/visit/statistics/plot-conrad-stats.sh

# plot all tracks
${MEANIE3D_HOME}/visit/visualize_results/visualize_conrad_tracks.sh .
