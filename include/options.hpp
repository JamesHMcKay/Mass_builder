#ifndef OPTIONS_H
#define OPTIONS_H
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include "utils.hpp"

using namespace std;


struct Options
{
public:
bool counter_terms = false;
bool verbose = false;
int loop_order = 2;
string particle;
string model;
string diagram;
string particle_safe;
string input_list="";

int run_mode = 0;
/* options are:  
0 not specified, throw error
1 calculate amplitude for specific diagram, if list provided
1b calculate amplitude for list of diagrams, if no list (-i list.txt) provided use models/<model>/diagrams.txt

choose between 1 and 1b depending on the input available

4 generate code for list of diagrams, if no list (-i list.txt) provided use models/<model>/diagrams_avail.txt
5 generate diagrams for all diagrams in list, if no list given generate all
6 generate diagram for specific 
*/



Options (){};
void print_options()
{
cout << "Choosen options are:" << endl;
cout << "counter terms = " << counter_terms << endl;
cout << "verbose = " << verbose << endl;
cout << "loop_order = " << loop_order << endl;
cout << "particle = " << particle << endl;
cout << "diagram = " << diagram << endl;
cout << "model = " << model << endl;
cout << "run mode = " << run_mode << endl;
cout << "input list = " << input_list << endl;
}
};

class User_input
{
public:

std::string param [100];
std::string name [100]; int n=0;
Options options;

User_input(int argc, char* argv[])
{
if (argc==1)
{
cout << "NO INPUT -- user input guide here" << endl;
}
else
{

for (int i = 0;i<argc;i++)
{
param[i] = argv[i];
}
n = argc;
}

}


bool find_string(string input)
{
bool result = false;
for (int i = 0;i<n;i++)
{
string s1 = param[i];
if (s1.find(input) != std::string::npos)
{
result = true;
}
}
return result;
}

bool find_and_read_string(string input,string &output)
{
bool result = false;
for (int i = 0;i<n;i++)
{
string s1 = param[i];
if (s1.find(input) != std::string::npos)
{
if ((i+1) >= n){cout << "flag to input " << output << " entered but no input specified after -- programme shutting down." << endl;}
else
{
 output = param[i+1];
 result = true;
}
}
}
return result;
}


void user_interface()
{

string io;

if (find_string("-c")){ options.counter_terms = true;}
if (find_string("-v")){ options.verbose = true;}
if (find_string("-a")){ options.run_mode = 1;}

if (find_string("-l"))
{
string input = "loop order";
if (find_and_read_string("-l",input))
{
if (input == "2"){ options.loop_order = 2;}
if (input == "1"){ options.loop_order = 1;}
else {cout <<"This loop order is not supported please enter 1 or 2"<<endl;}
}
}

if (find_string("-m"))
{
string input = "model name";
if (find_and_read_string("-m",input)){options.model = input;}
}

if (find_string("-p"))
{
string input = "particle name";
if (find_and_read_string("-p",input)){options.particle = input;}
}


if (find_string("-d"))
{
string input = "diagram number";
if (find_and_read_string("-d",input)){options.diagram = input;}
}


if (find_string("-i"))
{
string input = "a diagram list";
if (find_and_read_string("-i",input)){options.input_list = input;}
}



if (find_string("-g"))
{
string input = "generate code";
if (find_and_read_string("-g",input))
{
options.run_mode = 4;
options.model = input;
}
}

if (find_string("-gs"))
{
string input = "generate code";
if (find_and_read_string("-gs",input))
{
options.run_mode = 3;
options.model = input;
}
}

options.print_options();
}



};













#endif
