#ifndef SUPPLEMENTS_H
#define SUPPLEMENTS_H


#include "data.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>

#include <fstream>


#include "data.hpp"


namespace extra_TSIL_interface {


using namespace std;
typedef complex<double> dcomp;

class One_loop_derivatives
{
  private:
    
  public:
  Data data;
  One_loop_derivatives (){}
  One_loop_derivatives(Data _data)
  {
      data=_data;
   }
  
  
   double add_derivatives(Data &data);
   double add_derivatives_2(Data &data);
  
};
}

#endif