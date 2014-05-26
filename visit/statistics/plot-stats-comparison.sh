#!/usr/local/bin/gnuplot


# track length

#set terminal postscript eps enhanced color font 'Helvetica,12'
#set output "lengths-hist-comparison.eps"

set terminal svg enhanced
set output "lengths-hist-comparison.svg"
set style fill transparent solid 0.5
set title "Distribution of track length"
set xlabel "Track length in [#steps]"
set ylabel "Number of tracks"
set xtics 5
plot "lengths-hist.txt" with boxes title "M3D", "conrad_lengths-hist.txt" with boxes title "CONRAD"

# cluster size
set output "sizes-hist-comparison.eps"
set title "Distribution of cluster size"
set xlabel "log(cluster size in [#gridpoints])"
set ylabel "log(number of clusters)"
set logscale y
set logscale x
set xtics auto
plot "sizes-hist.txt" with boxes title "M3D", "conrad_sizes-hist.txt" with boxes title "CONRAD"
unset logscale y
unset logscale x

# speed
set output "speeds-hist-comparison.eps"
set title "Distribution of cluster speeds"
set xlabel "Cluster speed in [m/s]"
set ylabel "Number of clusters"
set xtics 2
plot "speeds-hist.txt" with boxes title "M3D", "conrad_speeds-hist.txt" with boxes title "CONRAD"

# directions
set output "directions-hist-comparison.eps"
set title "Distribution of tracking direction"
set xlabel "Cluster direction in [deg]"
set ylabel "Number of clusters"
set xtics 15
plot "directions-hist.txt" with boxes title "M3D", "conrad_directions-hist.txt" with boxes title "CONRAD"

