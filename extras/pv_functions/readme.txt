6th March 2016

This is simple program for calculating the MDM mass splitting using analytic expressions for the Passarino-Veltman
functions from hep-ph/9606211.

To compile:
> mkdir build
> cd build
> cmake ..
> make

There is one executable, main, run this and it will print out some example mass splitting calculations from
different methods, and it will plot one figure.  This figure will be found in Figures/Figures directory, and
the data used to generate it in Python is in the Figures/data directory.

Edit src/main.cpp to enable more figures or change the input values.  To edit the figures themselves edit
src/figures.cpp.