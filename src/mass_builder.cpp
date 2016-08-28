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

#define RUN_ALL

using namespace std;


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
else
{
run_calc_diagram(argc, argv);
}

}
