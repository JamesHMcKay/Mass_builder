#ifndef EW_TRIPLET_H
#define EW_TRIPLET_H

#include "data.hpp"


class EW_triplet_spectrum
{
private:
		
public:
  
  // Mass Builder data structure
  Data data;
   
  EW_triplet_spectrum(Data data) : data(data) {};
  
  EW_triplet_spectrum(){};
  
  // run FlexibleSUSY and update data structure with couplings and self energies
  bool compute_spectra_flexiblesusy();
  
  // set parameters according to method in Ibe et al. (2013)
  void compute_spectra_MB_1loop();
  void compute_spectra_MB_2loop();
  
  // compute self energies using TSIL
  void compute_tsil();
  void compute_tsil_iterative();
  
  // compute the iterative pole mass
  double iterative_ms_bar_mass(Data data, string particle);
  
  double get_deltam();
  double get_deltam_2loop();
  
  double get_charged_mass();
  
  double get_neutral_mass();
  
};

#endif
