#!/usr/bin/gnuplot
#
# For plotting a line chart.
#
# AUTHOR: Riccardo Battistini

reset

unset colorbox

# set background of the plot
set table 'thickbg.dat'
sp [0:5] [0:30] y
unset table
set grid front
unset grid
set palette defined (0 0.95 0.95 0.95, 1 0.95 0.95 0.95)

set key left maxrows 4

set xlabel "Time"
set ylabel "Price"
set tics font "Helvetica, 14" nomirror

set macros
background="[0:5] [0:30] 'thickbg.dat' w ima t ''"

set macros
dummy="NaN title ' ' lt -3"

plot @background
