#ifndef UTILS_H
#define UTILS_H


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

void print_math_body(ofstream &file,int loop_order,string particle_full,string diagram,string model,string cwd);

void print_product(ofstream &myfile,string name_1,string name_2,string SEn="SEn");

void print_A(ofstream &myfile, string elements,string SEn="SEn");

void print_B(ofstream &myfile, string elements,string SEn="SEn");

void print_V(ofstream &myfile, string elements,string SEn="SEn");

void print_T(ofstream &myfile, string elements,string SEn="SEn");

void print_J(ofstream &myfile, string elements,string SEn="SEn");

void print_K(ofstream &myfile, string elements,string SEn="SEn");

void print_F(ofstream &myfile, string elements,string SEn="SEn");

bool check_done();


}

#endif