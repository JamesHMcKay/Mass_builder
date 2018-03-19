#ifndef DECAYS_H
#define DECAYS_H


#include "pv.hpp"
#include "mdm.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>

#include <fstream>
using namespace std;
typedef complex<long double> dcomp;

class Decays
{
  private:
  Data data;
  MDM mdm;
  int method;
  long double deltam;

  public:
  Decays (){}  // defualt constructor
  Decays(Data _data,MDM _mdm)
  {
  data=_data;
  mdm=_mdm;
  method=data.method;
  
  //deltam=data.get_deltaM();
  } //constructor
  
   long double lifetime();
   long double lifetime_simple();
   long double lifetime2();
   long double lifetime_simple2();
   long double calc_lifetime(long double deltam);
  
   long double pion_channel(long double deltam);
  
   long double muon_channel(long double deltam);

   long double electron_channel(long double deltam);

   long double pion_decay();

  
};




#endif
