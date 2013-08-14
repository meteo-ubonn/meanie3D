#!/bin/bash
mkdir netcdf && mv *.nc netcdf
rm *_zeroshift_*.vtk
rm *_centers.vtk
rm *_modes.vtk
rm meanshift*.vtk
rm visitlog.py
mkdir vtk && mv *.vtk vtk
mkdir images && mv *.png images
mkdir log && mv *.log log
mkdir movies
convert -delay 50 -quality 100 images/source_*.png movies/source.mpeg
convert -delay 50 -quality 100 images/tracked_*.png movies/tracked.mpeg