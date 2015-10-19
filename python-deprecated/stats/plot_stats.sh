#!/usr/local/bin/gnuplot

# track length
set term postscript
set output "lengths-hist.eps"
set title "Distribution of track length"
set xlabel "Track length in [#steps]"
set ylabel "Number of tracks"
set xtics 5
plot "lengths-hist.txt" with boxes

# cluster size
set output "sizes-hist.eps"
set title "Distribution of cluster size"
set xlabel "log(cluster size in [#gridpoints])"
set ylabel "log(number of clusters)"
set logscale y
set logscale x
set xtics auto
plot "sizes-hist.txt" with boxes
unset logscale y
unset logscale x

# speed
set output "speeds-hist.eps"
set title "Distribution of cluster speeds"
set xlabel "Cluster speed in [m/s]"
set ylabel "Number of clusters"
set xtics 2
plot "speeds-hist.txt" with boxes

# directions
set output "directions-hist.eps"
set title "Distribution of tracking direction"
set xlabel "Cluster direction in [deg]"
set ylabel "Number of clusters"
set xtics 15
plot "directions-hist.txt" with boxes

