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

export scales="10 25 50 100 250"
echo "Running complete sets on scales ${scales}"

TRACKSTATS_PARAMS =  " --vtk-dimensions x,y"
TRACKSTATS_PARAMS += " --write-gnuplot-files"
TRACKSTATS_PARAMS += " --create-length-statistics --create-speed-statistics --create-direction-statistics --create-cluster-statistics"
TRACKSTATS_PARAMS += " -b raa"

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
    ${MEANIE3D_HOME}/visit/tracking/run_tracking-2D.sh "${1}" "${scale}"
    meanie3D-trackstats ${TRACKSTATS_PARAMS} -p netcdf
    rm -f nohup.out
    cd ..

done