#!/usr/bin/gnuplot
#
# For plotting figure 1.
#
# Last modified on 2021-04-09
# by Riccardo Battistini <riccardo.battistini(at)studio.unibo.it

reset

set loadpath './scripts'
load "histogram.gp"

set terminal svg size 900,400 rounded font "Open Sans, 13"
set xrange[0:8.95]
set yrange[0:105]

set key center top horizontal outside maxrows 4

set xlabel "Processi"
set ylabel "Wall clock time (s)"

# Change tics
set tics font "Open Sans, 11" nomirror out scale 0.75

#  Major tics
set xtics scale 0.01 out
set ytics 8 scale 0.01 in

# Minor tics
#set mxtics
#set mytics 2

plot "./datasets/mpi-test5-all.csv" u 2 ls 2 t "test5",\
"./datasets/mpi-test6-all.csv" u 2 ls 3 t "test6",\
"./datasets/mpi-test7-all.csv" u 2 ls 4 t "test7"
