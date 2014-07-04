#!/bin/bash

if [ "X$1" = "X" ]
then
    echo "run_tb4_clustering.sh <sourcedir>"
    echo "Creates a python script for cluster creation and tracking and runs it in Visit"
exit 0
fi

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
exit 0
fi

# clean up previous runs

if [ ! -d "images" ]
then
    mkdir images
else
    rm -rf images/*.png
fi

if [ ! -d "clusters" ]
then
    mkdir clusters
else
    rm -rf clusters/*.nc
fi

MEANIE3D_PARAMS="-d lat,lon --vtk-dimensions lon,lat -r 0.5,0.5,300 --weight-function=inverse --wwf-lower-threshold=0.001 -m 5 --write-weight-function"

FILES=`find $1 -name "*.nc"`
for file in $FILES
do
    fn=${file##*/}
    basename=${fn%.nc}

    echo "Processing file $file"

    for var in observation simulation
    do

        outfile="${basename}-${var}-clusters.nc"
        PARAMS="${MEANIE3D_PARAMS} -v $var --upper-thresholds $var=253.15 -f $file -o $outfile"
        echo "meanie3D-detect $PARAMS"
        meanie3D-detect $PARAMS

        echo "Visualizing results ..."
        echo "$MEANIE3D_HOME/visit/tb4/visualize_tb4_2D.sh $file $outfile $var"
        $MEANIE3D_HOME/visit/tb4/visualize_tb4_2D.sh $file $outfile $var

        exit

        mv ${basename}*.png images
    done

done

mv *.nc clusters
rm visitlog.py

