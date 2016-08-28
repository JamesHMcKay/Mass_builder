#!/bin/sh
# This should be run in root Mass_builder directory

#cd <system_dependent>/Mass_builder/

echo "Please enter path to the file mathematica script or use one of the default options below by entering the corresponding shortcut"
echo "  "
echo "osx    for OSX:/Applications/Mathematica.app/Contents/MacOS/MathematicaScript"
echo "osx64  for OSX Intel 64-bit machine:/Applications/Mathematica.app/Contents/MacOS/MathematicaScript64"
echo "math   for most other systems which already have kernal from /usr/local/bin/math"
echo "  "
read MATH_PATH
if [ "$MATH_PATH" == "osx" ]; then
  MATH_PATH="/Applications/Mathematica.app/Contents/MacOS/MathematicaScript"
fi

if [ "$MATH_PATH" == "osx64" ]; then
  MATH_PATH="/Applications/Mathematica.app/Contents/MacOS/MathematicaScript64"
fi

if [ "$MATH_PATH" == "math" ]; then
  MATH_PATH="/usr/local/bin/math"
fi

sed -i '' -e "s|.*MATH_PATH.*| /\*MATH_PATH \*/  file<< \"#!$MATH_PATH -script\"<<endl;|g" src/utils.cpp


echo "Please enter path to TSIL header tsil_cpp.h"
echo "for example /Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h"
read TSIL_PATH

if [ "$TSIL_PATH" == d ]; then
  TSIL_PATH=/Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h
fi


sed -i '' -e "s|.*TSIL_INCLUDE_PATH.*| /\*TSIL_INCLUDE_PATH \*/#include \"$TSIL_PATH\"  // Required TSIL header file|g" src/*.cpp

sed -i '' -e "s|.*TSIL_PATH.*| /\* TSIL_PATH \*/ std::string TSIL = \"$TSIL_PATH\";  |g" src/*.cpp
















cp src/self_energy_bak.cpp src/self_energy.cpp
mkdir build
mkdir build/output
mkdir build/generator

