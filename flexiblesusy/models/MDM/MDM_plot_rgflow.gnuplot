set terminal x11

set title "MDM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MDM_rgflow.dat'

plot for [i=2:34+1] filename using 1:(column(i)) title columnhead(i)
