#!/bin/sh
# This should be run in root Mass_builder directory

#cd <system_dependent>/Mass_builder/

echo "Please enter path to the folder containing MathematicaScript and MathKernel, use one of the default options below by entering the corresponding shortcut"
echo "  "
echo "osx    for OSX:/Applications/Mathematica.app/Contents/MacOS"
echo "math   for most other systems which already have kernal from /usr/local/bin/math"
read MATH_PATH
if [ "$MATH_PATH" == "osx" ]; then
  MATH_PATH="/Applications/Mathematica.app/Contents/MacOS"
fi
if [ "$MATH_PATH" == "math" ]; then
  MATH_PATH="/usr/local/bin/math"
fi

# add the Mathematica path into the required locations in the source files
sed -i '' -e "s|.*MATH_PATH.*|    /\*MATH_PATH \*/  file<< \"#!$MATH_PATH/MathematicaScript -script\"<<endl;|g" src/templates.cpp

sed -i '' -e "s|.*MATH_KERNEL_PATH.*|    /\*MATH_KERNEL_PATH \*/  WSTPflags << \"-linkname \" << \"$MATH_PATH/MathKernel\" << \" -mathlink\";|g" include/compute_amp.hpp

echo "  "
echo "Please enter path to directory containing the TSIL header tsil_cpp.h"
echo "for example /Users/<user_name>/Programs/tsil-1.3"
read TSIL_PATH

# following is default for my system, could use find / -name tsil_cpp.h to search but this is slow
if [ "$TSIL_PATH" == d ]; then
  TSIL_PATH=/Users/jamesmckay/Documents/Programs/tsil-1.3
fi

# add the TSIL path into the required places
sed -i '' -e "s|.*TSIL_INCLUDE_PATH.*|  /\*TSIL_INCLUDE_PATH \*/#include \"$TSIL_PATH/tsil_cpp.h\"  // Required TSIL header file|g" src/*.cpp

sed -i '' -e "s|.*TSIL_PATH.*|    /\* TSIL_PATH \*/ std::string TSIL = \"$TSIL_PATH/tsil_cpp.h\";  |g" src/*.cpp

sed -i '' -e "s|.*set(TSIL_HEADER_FILE.*|set(TSIL_HEADER_FILE $TSIL_PATH/tsil_cpp.h ) |g" CMakeLists.txt

sed -i '' -e "s|.*LINK_DIRECTORIES(.*|LINK_DIRECTORIES($TSIL_PATH) |g" CMakeLists.txt




# create bare self_energy.cpp so compile works before actual code generated
cp src/self_energy_bak.cpp src/self_energy.cpp
cp include/self_energy_backup.hpp include/self_energy.hpp
cp include/data_bak.hpp include/data.hpp
#rm src/amp_*.cpp
# create working directories and build directory
mkdir build
mkdir output
mkdir models/EW_triplet/output
mkdir models/Scalar/output

