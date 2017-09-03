set terminal x11

set title "EW_triplet particle spectrum"
set ylabel "mass / GeV"
unset key
unset bars

if (!exists("filename")) filename='EW_triplet_spectrum.dat'

plot filename using 1:2:(0.4):xtic(3) with xerrorbars pointtype 0 linecolor rgb "black"
