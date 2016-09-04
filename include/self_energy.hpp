#ifndef SELF_ENERGY_H
#define SELF_ENERGY_H


#include "data.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>

#include <fstream>




using namespace std;

typedef complex<double> dcomp;

class Self_energy
{
  private:
  
  
  public:
  Data data;
  Self_energy (){}  // defualt constructor
  void run_tsil(Data &data);
  
  void init_tsil(Data data);
  
  
   
};
//}

#endif