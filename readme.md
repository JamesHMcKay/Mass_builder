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

This program requires Mathematica, FeynCalc, FeynArts (patched for use with FeynCalc) and TARCER to be installed.  For instructions on how to install these please see the user manual.  To install Mass Builder following the instructions below.

If using a Linux machine then first run the configure script
```
./scripts/config_linux.sh or 
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

