#!/usr/local/bin/gnuplot

# NOTE: this script expects the conrad_ - versions of the histogram
# files to be present in the directory ABOVE the one you're plotting
# in. This is due to the fact that Meanie3D runs on numerous scales
# and the CONRAD stats would have to be copied for each scale. It's
# easier this way.

set terminal png
set style fill transparent pattern 4 border 2

# track length
set title "Distribution of track length"
set xlabel "Track length in [#steps]"
set ylabel "Number of tracks"
set xtics 5
set xrange [2:50]
set output "lengths-hist-comparison.png"
plot "lengths-hist.txt" with boxes title "M3D", "../conrad_lengths-hist.txt" with boxes title "CONRAD"
set autoscale x

# cluster size
set title "Distribution of cluster size"
set xlabel "log(cluster size in [#gridpoints])"
set ylabel "number of clusters"
set logscale x
set xrange [10:10000]
set xtics auto
set output "sizes-hist-comparison.png"
plot "sizes-hist.txt" with boxes title "M3D", "../conrad_sizes-hist.txt" with boxes title "CONRAD"
unset logscale x
set autoscale x

# speed
set title "Distribution of cluster speeds"
set xlabel "Cluster speed in [m/s]"
set ylabel "Number of clusters"
set xtics 5
set output "speeds-hist-comparison.png"
plot "speeds-hist.txt" with boxes title "M3D", "../conrad_speeds-hist.txt" with boxes title "CONRAD"

# directions
set title "Distribution of tracking direction"
set xlabel "Cluster direction in [deg]"
set ylabel "Number of clusters"
set xtics 30
set xrange [15:360]
set output "directions-hist-comparison.png"
plot "directions-hist.txt" with boxes title "M3D", "../conrad_directions-hist.txt" with boxes title "CONRAD"

