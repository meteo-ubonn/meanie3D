#!/usr/local/bin/gnuplot

set terminal png
set style fill transparent pattern 4 border 2

# track length
set title "Distribution of track length"
set xlabel "Track length in [#steps]"
set ylabel "log(#tracks)"
set xtics 5
set xrange [2:50]
set logscale y
set output "lengths-hist-comparison.png"
plot \
"scale10/lengths-hist.txt" with linespoints title "t=10km",\
"scale25/lengths-hist.txt" with linespoints title "t=25km"\

unset logscale y
set autoscale x

# cluster size
set title "Distribution of cluster size"
set xlabel "log(cluster size in [#gridpoints])"
set ylabel "number of clusters"
set logscale x
set xrange [10:10000]
set xtics auto
set output "sizes-hist-comparison.png"
plot \
"scale10/sizes-hist.txt" with linespoints title "t=10km", \
"scale25/sizes-hist.txt" with linespoints title "t=25km" \

unset logscale x
set autoscale x

# speed
set title "Distribution of cluster speeds"
set xlabel "Cluster speed in [m/s]"
set ylabel "Number of clusters"
set xtics 5
set output "speeds-hist-comparison.png"
plot \
"scale10/speeds-hist.txt" with linespoints title "t=10km", \
"scale25/speeds-hist.txt" with linespoints title "t=25km" \

# directions
set title "Distribution of tracking direction"
set xlabel "Cluster direction in [deg]"
set ylabel "Number of clusters"
set xtics 30
set xrange [15:360]
set output "directions-hist-comparison.png"
plot \
"scale10/directions-hist.txt" with linespoints title "t=10km", \
"scale25/directions-hist.txt" with linespoints title "t=25km" \



