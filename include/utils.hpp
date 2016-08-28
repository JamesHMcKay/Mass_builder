#ifndef UTILS_H
#define UTILS_H

// useful functions that are super specific to any area of the code


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


using namespace std;

namespace utils
{

void get_data(vector<std::string> &A,int &n,const char *filename);

void get_data(vector<std::string> &A, vector<std::string> &B,int &n,const char *filename);

void print_math_header(ofstream &myfile);





}

#endif