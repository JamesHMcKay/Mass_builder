#ifndef SELF_ENERGY_H
#define SELF_ENERGY_H

#include "data.hpp"
#include <complex>


using namespace std;

typedef complex<double> dcomp;

class Self_energy
{
  private:
  
  public:
  Data data;
  Self_energy (){}
  void run_tsil(Data &data);
  
  void init_tsil(Data data);
   
};

#endif