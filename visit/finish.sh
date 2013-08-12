#!/bin/bash
mkdir netcdf && mv *.nc netcdf
mkdir vtk && mv *.vtk vtk
mkdir images && mv *.png images
mkdir log && mv *.log log
mkdir movies
convert -delay 50 -quality 100 images/source_*.png movies/source.mpeg
convert -delay 50 -quality 100 images/tracked_*.png movies/tracked.mpeg