#ifndef CALC_AMPLITUDES_H
#define CALC_AMPLITUDES_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"

using namespace utils;

class Calc_amplitudes
{
public:

Calc_amplitudes(){}

bool calc_diagram(Options options);
void generate_figures(Options options);
};

#endif