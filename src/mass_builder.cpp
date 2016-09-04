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
using namespace utils;

Options options;



// main routine to manage a diagram by diagram procedure
void run_mass_builder_mode_1a(Options options)
{
Calc_amplitudes ca;

std::string particles [1000];
std::string types [1000];
std::string diagrams [1000]; int i=0;
std::string model = options.model;
const char *file_diagrams;
if (options.input_list==""){
const char *ext = ".txt";
const char* file_diagrams_tmp = "models/";
string c_file_diagrams = file_diagrams_tmp + model + "/diagrams" + ext;
file_diagrams = c_file_diagrams.c_str();
}
else
{
file_diagrams = options.input_list.c_str();
}


std::ifstream input(file_diagrams);
std::string line;
while(getline(input, line))
{
  if (!line.length() || line[0] == '#')
     continue;
  std::istringstream iss(line);
  iss>> particles[i] >> diagrams[i] >> types[i];
  i=i+1;
}
// run over all entries in the input file
for (int k=0;k<i;k++)
{
options.particle = particles[k];
options.diagram = diagrams[k];
options.set_type(types[k]);
options.model = model;
ca.calc_diagram(options);
}
}





void run_mass_builder_mode_1b(Options options)
{
Calc_amplitudes ca;
ca.calc_diagram(options);
}



void evaluate(Options options)
{

Self_energy se;
Data data(options);
se.run_tsil(data);


}



int main(int argc, char *argv[])
{

User_input user(argc,argv);

user.user_interface();

Options options = user.options;

// read options and work through possibilities for each run mode and check requirements are meant
if (options.model=="" && (options.run_mode != 6)){ cout << "no model specified" << endl; return 0;}


if (options.run_mode == 1)
{
// going to calculate amplitudes
if ((options.particle == "") || (options.diagram == ""))
{
run_mass_builder_mode_1a(options);
}
else
{
if ((options.input_list==""))
{
if ((options.particle == "") || (options.diagram == "")) { cout << "no valid input selected" << endl;}
else { run_mass_builder_mode_1b(options);}
}
}
}


if (options.run_mode == 4 )
{
if (options.input_list == "")
{
options.input_list = "models/"+ options.model + "/output/avail_diagrams_.txt";
}
if (options.model == "") { cout << "please specify a model" << endl; return 0;}
Generate_code::generate_code(options);

}


if (options.run_mode == 5 )
{
if (options.model == "" || options.particle == "") { cout << "please specify a model and particle, at least one is missing" << endl; return 0;}
Calc_amplitudes ca;
ca.generate_figures(options);

}


if (options.run_mode == 6 )
{
if (options.input_list == "" ) { cout << "missing input data" << endl; return 0;}

evaluate(options);

}










/*
if (argc == 1)
{
user_input_guide();
}
else
{
string option = argv[1];
if (option=="-g")
{
Generate_code::main_function(argc,argv,options);
}
else if (option=="-gc")
{
options.counter_terms = true;
Generate_code::main_function(argc,argv,options);
}
else if (option == "-d")
{
Calc_amplitudes ca;
ca.generate_figures(argc,argv,options);
}
else if (option == "-dc")
{
options.counter_terms = true;
Calc_amplitudes ca;
ca.generate_figures(argc,argv,options);
}
else if (option == "-e")
{
evaluate(argc, argv);
}
else if (option == "-a")
{
run_calc_diagram(argc, argv);
}
else if (option == "-ac")
{
options.counter_terms = true;
run_calc_diagram(argc, argv);
}
else
{

user_input_guide();


}
}*/


}
