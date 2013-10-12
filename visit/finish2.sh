#!/bin/bash
mkdir netcdf && mv *.nc netcdf
mkdir log && mv *.log log
rm -f *.vt*
rm visitlog.py
meanie3D-trackplot --basename raa01-rx_10000-130620 --sourcepath netcdf --create-cumulated-size-statistics 0