Description:

This program is designed to completely automate the process of two-loop self energy calculation taking an input of a FeynArts model file all the way through to a numerical result for the self energy.

We make use of the following tools to complete this process:

FeynArts (hep-ph/0012260) is used to generate the two-loop amplitudes
FeynCalc (ArXiv:1601.01167) is used to reduce the amplitudes
TARCER (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
TSIL (hep-ph/0501132) is used to evaluate the basis integrals

Required programmes:

To use this tool you need to have:

Mathematica
FeynCalc including a FeynArts installation (ships with by default)  # add more details here
TSIL


Instructions:

Quick example:

cd to Mass_builder directory

./scripts/config.sh
cd build
cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ -DCMAKE_C_COMPILER=/usr/local/bin/gcc ..
make
cd ..
./mass_builder -f Scalar      #  this step will calculate all the amplitudes specified in models/Scalar/diagrams.txt
./mass_builder -g Scalar      #  generate the TSIL interface code in src/self_energy.cpp
cd build
make					#  must rebuild to compile TSIL interface
cd ..
./mass_builder -e input.txt  #  evaluate self energy using the model parameters specified in input.txt