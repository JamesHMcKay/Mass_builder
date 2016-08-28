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


namespace supplementary_code {


using namespace std;
typedef complex<double> dcomp;

class Supplements
{
  private:
    
  public:
  Data data;
  Supplements (){}  // defualt constructor
  Supplements(Data _data)
  {
      data=_data;
   } //constructor
  
  
   void add_derivatives(/*Data &data*/);
  
  
  
   
};
}

#endif