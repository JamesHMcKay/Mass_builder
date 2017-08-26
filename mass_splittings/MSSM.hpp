#ifndef MSSM_H
#define MSSM_H

#include "data.hpp"


class MSSM
{
private:
  
public:
  
  MSSM(){};
  
  // run FlexibleSUSY and update data structure with spectrum
  void set_flexiblesusy_spectra(Data &data);
  
  // set parameters according to method in Ibe et al. (2013)
  void set_SM_parameters_1loop(Data &data);
  void set_SM_parameters_2loop(Data &data);
  
  // temporary for now, get difference of 1-loop derivatives
  double add_1loop_der(Data &data);
  
};

#endif
