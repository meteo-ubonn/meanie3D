#!/bin/bash

# Check some pre-requisites

if [ $# -lt 1 ]
then
    echo "create_all_oase_loops_2D <netcdf directory>"
    echo "Creates loops of all oase composite variables and a composite loop"
    exit 0
fi

# scripts reachable?

if [ "X${MEANIE3D_HOME}" = "X" ]
then
    echo "Please set environment variable MEANIE3D_HOME"
    exit 0
fi

# netcdf directory exists?

if [ ! -d "$1" ]; then
    echo "Can't find input directory $1"
    exit 0
fi

# Visualize variables individually

variables=(cband_radolan_rx linet_oase_tl msevi_l15_ir_108 msevi_l15_vis006 msevi_l2_cmsaf_cot msevi_l2_cmsaf_cph msevi_l2_cmsaf_cwp msevi_l2_cmsaf_reff msevi_l2_nwcsaf_cma msevi_l2_nwcsaf_crr msevi_l2_nwcsaf_ct msevi_l2_nwcsaf_cth)

for i in "${variables[@]}"
do
    ${MEANIE3D_HOME}/visit/loops/create_oase_loop_2D.sh $1 $i
done

# Run the composite 

${MEANIE3D_HOME}/visit/loops/create_composite_loop_2D.sh $1