#!/bin/bash

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
exit 0
fi

export dirs=`ls -d scale*`
echo "Running tracking statistics on dirs: $dirs"

for dir in $dirs; do
    if [ -d $dir/netcdf ];then
        cd $dir
        meanie3D-trackstats -t -g -t -1 -2 -3 -4 --vtk-dimensions=x,y -s netcdf
        ${MEANIE3D_HOME}/scripts/stats/plot_stats.sh
        ${MEANIE3D_HOME}/scripts/visit/tracks/visualize_radolan_tracks.sh .
        cd ..
    fi
done