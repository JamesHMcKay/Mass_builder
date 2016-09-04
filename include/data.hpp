#ifndef DATA_H
#define DATA_H
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "options.hpp"
using namespace std;
struct Data
{
public:
double lambda;
double dc;
double dg;
double Ms;
double SE_S1;
double P, Q;
  Data (){};
Data(Options options) {
double param [99];
std::string name [99]; int i=0;
std::ifstream input(options.input_list);
std::string line;
while(getline(input, line)) {
if (!line.length() || line[0] == '#')
 continue;
 std::istringstream iss(line);
 iss>> name[i] >> param[i];
     i=i+1;
   }
  for (int n=0;n<i+1;n++)
  {
  if (name[n]=="lambda")
  {
  lambda=param[n];
  }
  if (name[n]=="dc")
  {
  dc=param[n];
  }
  if (name[n]=="dg")
  {
  dg=param[n];
  }
  if (name[n]=="Ms")
  {  Ms=param[n];  }
  if (name[n]=="Q")
  {  Q =param[n];  }
  if (name[n]=="P")
  {  P =param[n];  }
}
}
};
#endif
