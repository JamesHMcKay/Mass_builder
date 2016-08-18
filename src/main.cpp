/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 18/09/16
James McKay

--- main.cpp ---

Reads the user input file and initialises the Data structure
then passes this onto the run_tsil function.

*/

#include "self_energy.hpp"
#include "data.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>

using namespace std;



int main(int argc, char* argv[])
{

Data data(argc,argv);

Self_energy se(data);


se.run_tsil(data);


cout << "SE_1 - SE_2 = " << data.SE_1 - data.SE_2 << endl;

}