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
            -b raa -p netcdf
        ${MEANIE3D_HOME}/visit/statistics/plot-stats.sh
        cd ..
    fi
done