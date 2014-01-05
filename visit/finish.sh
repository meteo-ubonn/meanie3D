#!/bin/bash
mkdir netcdf && mv *.nc netcdf
rm -f *_zeroshift_1*.vtk
rm -f *_zeroshift_2*.vtk
rm -f *_zeroshift_3*.vtk
rm -f *_zeroshift_4*.vtk
rm -f *_zeroshift_*.vtk
rm -f *_centers.vtk
rm -f *_modes.vtk
rm -f meanshift*.vtk
rm visitlog.py
mkdir vtk 
mv *.vtk vtk
mv *.vtr vtk
mkdir images && mv *.png images
mkdir log && mv *.log log
mkdir movies
convert -limit memory 4GB -delay 50 -quality 100 images/source_*.png movies/source.mpeg
convert -limit memory 4GB -delay 50 -quality 100 images/tracked_*.png movies/tracked.mpeg
meanie3D-trackstats --basename raa01-rx_10000-130620 --sourcepath netcdf --create-cumulated-size-statistics --create-direction-statistics