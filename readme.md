Mass Builder
======

This program is designed to completely automate the process of two-loop self energy calculation taking an input of a FeynArts model file all the way through to a numerical result for the self energy.

We make use of the following tools to complete this process:

- **FeynArts:** (hep-ph/0012260) is used to generate the two-loop amplitudes
- **FeynCalc:**  (ArXiv:1601.01167) is used to reduce the amplitudes
- **TARCER:**  (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
- **TSIL:**  (hep-ph/0501132) is used to evaluate the basis integrals


Documentation
--
This is an interface tool making use of the existing Mathematica and C packages to compute a numerical self energy.  The process is separated into three main steps, although many other features are available.  1) Amplitude calculation — this is where we run FeynCalc and and decompose the amplitude into a list of basis integrals and corresponding coefficients. 2) Code Generation — in this step we used the stored output from step one to generate C++ code which can interface to the TSIL libraries.  3) Code evaluation — this is the numerical evaluation of the self energy using the TSIL libraries.



Installation
--

This program requires Mathematica, FeynCalc, FeynArts (patched for use with FeynCalc) and TARCER to be installed.  For instructions on how to install these please see the user manual.  

The easiest way to install FeynCalc, FeynArts and Tarcer is via the automated installation method.  Open a Mathematica notebook or kernel session and enter
```
Import["https://raw.githubusercontent.com/FeynCalc/feyncalc/master/install.m"]
InstallFeynCalc[]
```
when requested to install the latest version of FeynArts say yes, as this will automatically patch the FeynArts installation.  If you do not follow this method then it is not possible to run FeynArts and FeynCalc in the same session (as we need to do) as many function names are identical between the packages, so to avoid name shadowing follow the recommend method.

In the same notebook session run the following command to generate the Tarcer files
```
GenerateTarcerMX
```
All packages within Mathematica are now set up.

To install the Two-loop Self-energy Integral Library (TSIL) downloaded the source from http://www.niu.edu/spmartin/TSIL/. It may installed anywhere (as Mass Builder will request the path at configuration).

To install Mass Builder follow the instructions below.

If using a Linux machine then first run the configure script
```
./scripts/config_linux.sh
```
or for OS X run
```
./scripts/config.sh
```
this will prompt for the location of the MathematicaScript, and the folder containing TSIL.

Next to compile the code run

```
cd build
cmake ..
make
```

Quick start guide — basic example
--

To test the program is functioning correctly run the following example.


```
./mass_builder -a -m Scalar
./mass_builder -g -m Scalar
cd build
make
cd ..
./mass_builder -e -i models/Scalar/input.txt
```

this will return the result for the self energy of the simple scalar field theory example

```
One loop self energy of particle S1 = -0.0316688
Two loop self energy of particle S1 = 2.91938e-05
```
see the user manual for more details of this model and other usage examples.

Advanced example
--
In this example we compute the one-loop self energies for an electroweak multiplet, in the context of the MSSM.  Following the commands below to generate a figure showing mass splittings as a function of the electroweak multiplet tree-level mass.

```
mkdir models/MSSM/output
./mass_builder -a -m MSSM -i models/MSSM/lists/example_1.txt
./mass_builder -g -m MSSM -i models/MSSM/lists/example_1.txt
cd build
cmake .
make MSSM
```
This will build a new executable that demonstrates how one may call Mass Builder from external functions.  The source file is located in examples/MSSM.cpp.  Once this has build run the following to generate a figure.

```
./MSSM -i models/MSSM/input.txt
python examples/plot_example.py
```
This will place a figure mass_splittings_MSSM.eps in the root directory.

