#ifndef MSSM_H
#define MSSM_H

#include "data.hpp"


class MSSM
{
private:
		
public:
  
  // Mass Builder data structure
  Data data;
   
  MSSM(Data data) : data(data) {};
  
  // run FlexibleSUSY and update data structure with couplings and self energies
  void compute_spectra_flexiblesusy();
  
  // set parameters according to method in Ibe et al. (2013)
  void compute_spectra_MB_1loop();
  void compute_spectra_MB_2loop();
  
  // compute self energies using TSIL
  void compute_tsil();
  void compute_tsil_iterative();
  
  // temporary for now, get difference of 1-loop derivatives
  double add_1loop_der(Data &data);
  
};

#endif
