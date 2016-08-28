Description:

This program is designed to completely automate the process of two-loop self energy calculation taking an input of a FeynArts model file all the way through to a numerical result for the self energy.

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