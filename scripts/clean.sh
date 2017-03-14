#!/bin/sh
# This should be run in root Mass_builder directory

# create bare self_energy.cpp so compile works before actual code generated
cp src/self_energy_bak.cpp src/self_energy.cpp
cp include/self_energy_backup.hpp include/self_energy.hpp
cp include/data_bak.hpp include/data.hpp
rm src/amp_*.cpp