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

#include "calc_amplitudes.hpp"
#include "generate_code.hpp"

#include "self_energy.hpp"
#include "supplements.hpp"

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

#define RUN_ALL

using namespace std;
using namespace supplementary_code;

// main routine to manage a diagram by diagram procedure
void run_calc_diagram(int argc, char *argv[])
{
Calc_amplitudes ca;

string diagram = "";
string particle = "";
string model = "";
if (argc==1)
{
cout << "Please enter diagram number to calcuate and particle name, options are \"chi0\" and \"chi1\"  " << endl;
cout << "or alternatively enter the name to a list as \"-f list.txt\" " << endl;
}
else
{

string option = argv[1];
if (option=="-f")
{
  std::string particles [100];
  std::string diagrams [100]; int i=0;
  std::ifstream input(argv[2]);
  std::string line;
  while(getline(input, line))
  {
        if (!line.length() || line[0] == '#')
           continue;
        std::istringstream iss(line);
        iss>> particles[i] >> diagrams[i];
        i=i+1;
  }
// run over all entries in the input file
  for (int k=0;k<i;k++)
  {
  
  ca.calc_diagram(diagrams[k],particles[k],argv[3]);
  }

}
else
{
model = argv[3];
cout << "using model = " << model << endl;
diagram = argv[2];
particle = argv[1];
ca.calc_diagram(diagram,particle,model);
}


}
}

void evaluate(int argc, char *argv[])
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



int main(int argc, char *argv[])
{

string option = argv[1];
if (option=="-g")
{
Generate_code::main_function(argc,argv);
}
else if (option == "-a")
{
Calc_amplitudes ca;
ca.generate_figures(argc,argv);
}
else if (option == "-e")
{
// evaluate self energy
evaluate(argc, argv);
}
else
{
run_calc_diagram(argc, argv);
}

}
