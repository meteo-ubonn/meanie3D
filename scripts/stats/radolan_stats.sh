#!/bin/bash

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

meanie3D-trackstats  meanie3D-trackstats -g -t -1 -2 -3 -4 --vtk-dimensions=x,y -s netcdf
${MEANIE3D_HOME}/scripts/stats/plot_stats.sh
${MEANIE3D_HOME}/scripts/visit/tracks/visualize_radolan_tracks.sh .