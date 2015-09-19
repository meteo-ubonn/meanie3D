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
        ${MEANIE3D_HOME}/scripts/stats/plot-stats-comparison.sh
        cd ..
    fi
done