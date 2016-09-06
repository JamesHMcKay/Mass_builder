#ifndef OPTIONS_H
#define OPTIONS_H
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>



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
5 draw diagrams for all diagrams in list, if no list given generate all

6 evaluate self energy using input in model folder
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


void set_type(std::string type);


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
call_user_guide();
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

void call_user_guide();

void user_interface();

};






#endif
