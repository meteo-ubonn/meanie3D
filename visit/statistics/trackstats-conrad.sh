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
    -f $1

${MEANIE3D_HOME}/visit/statistics/plot-conrad-stats.sh
