/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 18/09/16
James McKay

FeynArts (hep-ph/0012260) is used to generate the two-loop amplitudes
FeynCalc (ArXiv:1601.01167) is used to reduce the amplitudes
TARCER (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
TSIL (hep-ph/0501132) is used to evaluate the basis integrals


--- main.cpp ---

Reads the user input file and initialises the Data structure
then passes this onto the run_tsil function.

*/



#include "data.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>

#include "self_energy.hpp"
#include "supplements.hpp"

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

using namespace std;

using namespace supplementary_code;

int main(int argc, char* argv[])
{

Data data(argc,argv);

Self_energy se(data);

se.run_tsil(data);

/*

//se.run_tsil(data);
cout << "SE_1 - SE_2 = " << data.SE_1 - data.SE_2 << endl;

// FeynArts assigns a loop factor of -1/(16.0L*PI^4)^2 for two loop amplitudes
// but we should be using 1/(16.0L*PI^2)^2,
// so we rescale by PI^4 where appropriate
data.SE_1 = data.SE_1*pow(PI,4);
data.SE_2 = data.SE_2*pow(PI,4);

Supplements supp;
supp.add_derivatives(data);

cout << "SE_1 - SE_2 = " << data.SE_1 - data.SE_2 << endl;
*/

}
