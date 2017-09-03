#ifndef MSSM_SPECTRUM_H
#define MSSM_SPECTRUM_H

#include "data.hpp"


class MSSM_spectrum
{
private:
		
public:
  
  // Mass Builder data structure
  Data data;
   
  MSSM_spectrum(Data data) : data(data) {};
  
  MSSM_spectrum(){};
  
  // run FlexibleSUSY and update data structure with couplings and self energies
  bool compute_spectra_flexiblesusy();
  
  // set parameters according to method in Ibe et al. (2013)
  void compute_spectra_MB_1loop();
  void compute_spectra_MB_2loop();
  
  // compute self energies using TSIL
  void compute_tsil();
  void compute_tsil_iterative();
  
  // temporary for now, get difference of 1-loop derivatives
  double add_1loop_der(Data &data);
  
  // compute the iterative pole mass
  double iterative_ms_bar_mass(Data data, string particle);  
  
  double get_deltam();
  double get_deltam_2loop();
  
  double get_charged_mass();
  
  double get_neutral_mass();
  
};

#endif