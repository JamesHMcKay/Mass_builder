#!/bin/sh
# This should be run in root Mass_builder directory

if [ -f "src/self_energy.cpp" ]; then
  rm src/amp_*.cpp
fi

# create bare self_energy.cpp so compile works before actual code generated
cp src/self_energy_bak.cpp src/self_energy.cpp
cp include/self_energy_backup.hpp include/self_energy.hpp
cp include/data_bak.hpp include/data.hpp

# create working directories and build directory

if [ ! -d "build" ]; then
  mkdir build
fi
if [ ! -d "output" ]; then
  mkdir output
fi
if [ ! -d "models/EW_triplet/output" ]; then
  mkdir models/EW_triplet/output
fi
if [ ! -d "models/Scalar/output" ]; then
  mkdir models/Scalar/output
fi
if [ ! -d "models/MSSM/output" ]; then
  mkdir models/MSSM/output
fi
if [ ! -d "models/VDM/output" ]; then
  mkdir models/VDM/output
fi