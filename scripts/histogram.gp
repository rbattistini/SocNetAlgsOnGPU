#!/usr/bin/gnuplot
#
# For plotting a histogram.
#
# Last modified on 2021-04-09
# by Riccardo Battistini <riccardo.battistini(at)studio.unibo.it

#set terminal svg size 900,400 rounded font "Open Sans, 13"
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
set style fill solid 1.0
set datafile separator ","

unset key        # toggle legend

set xlabel "Processi"
set ylabel "Wall clock time (s)"

#  Major tics
set ytics
set xtics

# Minor tics
#set mxtics
#set mytics 2

# Change tics
set tics nomirror out scale 0.75

# Line styles for grid (81 & 102)
set style line 81 lt 0 lc rgb "#1F2430" lw 0.5
set style line 102 lc rgb'#808080' lt 0 lw 1

# Draw the grid lines for both the major and minor tics
set grid ytics back ls 81
set grid mytics back ls 81

# Line style for axes (80 & 101)
# set style line 80 lt 8 lc rgb "#808080"
set style line 101 lc rgb '#1F2430' lt 1 lw 1.3

set border 1 ls 101

# Line styles for data (1..6)
set style line 1 lc rgb "#58841D" # Vida Loca
set style line 2 lc rgb "#1F6EBB" # San Marino
set style line 3 lc rgb "#4E4E4E" # Dark grey
set style line 4 lc rgb "#DB7E54" # Copperfield
set style line 5 lc rgb "#914423" # Cumin
set style line 6 lc rgb "#73AF26" # Christi
