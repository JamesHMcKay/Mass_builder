Mass Builder
======

This program is designed to build, up from the level of a FeynArts model file, a C++ computer code to evaluate renormalised masses. This is achieved by generating the necessary Mathematica and C++ scripts to interface with the existing tools, along with sophisticated intermediary sorting.

I make use of the following tools to complete this process:

- **FeynArts:** (hep-ph/0012260) is used to generate the two-loop amplitudes
- **FeynCalc:**  (arXiv:1601.01167) is used to reduce the amplitudes
- **TARCER:**  (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
- **TSIL:**  (hep-ph/0501132) is used to evaluate the basis integrals
- **FIRE5:**  (arXiv:1408.2372) to reduce basis integrals
- **FeynHelpers** (arXiv:1611.06793) to communicate with FIRE5

I also make use of FindMathematica.cmake (Sascha Kratky 2010-2016 ) for locating the Mathematica installation within the cmake system.

Documentation
--
This is an interface tool making use of the existing Mathematica and C packages to compute a numerical self energy.  The process is separated into three main steps, although many other features are available.  1) Amplitude calculation — this is where we run FeynCalc and and decompose the amplitude into a list of basis integrals and corresponding coefficients. 2) Code Generation — in this step we used the stored output from step one to generate C++ code which can interface to the TSIL libraries.  3) Code evaluation — this is the numerical evaluation of the self energy using the TSIL libraries.


Installation
--

This program requires Mathematica, FeynCalc, FeynArts (patched for use with FeynCalc) and TARCER to be installed.

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

To install the Two-loop Self-energy Integral Library (TSIL) download the source from http://www.niu.edu/spmartin/TSIL/. It may installed anywhere (as Mass Builder will request the path at configuration).

To install Mass Builder enter the following commands

```
mkdir build
mkdir output
cd build
cmake -DTSIL_PATH=/path/to/tsil-1.4/ ..
make -jn
```
where you must specify the location of the TSIL directory as a flag to the cmake call , and n corresponds to the number of processes you have available.

Quick start guide — basic example
--

To test the program is functioning correctly run the following example.


```
mkdir models/Scalar/output
./mass_builder -a -m Scalar
./mass_builder -g -m Scalar
cd build
cmake ..
make
cd ..
./mass_builder -e -i models/Scalar/input.txt
```

this will return the result for the self energy of the simple scalar field theory example

```
One loop self energy of particle S1 = -0.0316688
Two loop self energy of particle S1 = 2.91938e-05
```

Before generating further code it is important to run the scripts/clean.sh again, as this will clean the previous generated files (such that cmake won't find these and try and compile them with the new, potentially inconsistent, code).

Supported models
--
I support five models with Mass Builder, although implementing new models is straight forward.  The available models are

- **Scalar:** a simple scalar field theory with a cubic and quartic interaction
- **MDM:** an electroweak quintuplet model
- **MSSM:** a modified version of the MSSM FeynArts model file shipped with FeynArts version 3.9
- **VDM:** a vector multiplet extension of the SM
- **QED:** a simple massless QED model

Further information
--

Please see the documentation in Mass_builder.pdf for detailed information on the algorithm structure and the many features available in code.  Some of these features are:

- convenient printing of FeynArts diagrams
- automatic computation of tree-level counter-term couplings
- generation of LaTeX ready list of Feynman rules for requested vertices
- use of parallel computing for many diagram computations
- mixed particle interactions such as Z -> W at one-loop
- optimisation of TSIL interface code by reducing function calls





