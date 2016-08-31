#ifndef CALC_AMPLITUDES_H
#define CALC_AMPLITUDES_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <ctime>

#include <unistd.h>
#include "utils.hpp"
#include "bases.hpp"

class Calc_amplitudes
{
public:

Calc_amplitudes(){}

bool calc_diagram(std::string diagram,std::string particle,std::string model);
void generate_figures(int argc, char *argv[]);
};

#endif