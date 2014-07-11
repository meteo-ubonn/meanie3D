#!/bin/bash

if [ "X$1" = "X" ]; then
    echo "multiscale_run_RX-2D.sh <path to directory containing composite files>"
    echo "Runs the clustering for a whole set of scales"
    exit 0
fi

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

export scales="5 10 25 50"
echo "Running complete sets on scales ${scales}"

#
# Iterate over scales
#
for scale in $scales; do

    export dest="scale${scale}"

    if [ ! -d $dest ]; then
        echo "Creating directory ${dest}"
        mkdir ${dest}
    else
        echo "Directory exists. Removing content"
        rm -rf ${dest}/*
    fi

    cd ${dest}
    ${MEANIE3D_HOME}/visit/tracking/run_oase_tracking-3D.sh "$1" "$scale"
    rm -f nohup.out
    cd ..

done

#
# Collate stats
#
${MEANIE3D_HOME}/visit/statistics/trackstats.sh

#
# Visualize
#
for scale in $scales; do

done
