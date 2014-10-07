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
        meanie3D-trackstats --write-track-dictionary --write-gnuplot-files \
            --create-length-statistics --create-speed-statistics \
            --create-direction-statistics --create-cluster-statistics \
            --write-center-tracks-as-vtk --vtk-dimensions=x y \
            -b raa -p netcdf
        ${MEANIE3D_HOME}/scripts/stats/plot-stats.sh
        ${MEANIE3D_HOME}/scripts/visit/tracks/visualize_radolan_tracks.sh .
        cd ..
    fi
done