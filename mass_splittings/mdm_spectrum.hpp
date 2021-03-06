#ifndef MDM_SPECTRUM_H
#define MDM_SPECTRUM_H

#include "data.hpp"


class MDM_spectrum
{
private:
		
public:
  
  // Mass Builder data structure
  Data data;
  double MChi_pole;
   
  MDM_spectrum(Data data) : data(data) {};
  
  MDM_spectrum(){};
  
  // run FlexibleSUSY and update data structure with couplings and self energies
  bool compute_spectra_flexiblesusy(int loop_order = 1, bool mass_ql_zero = false);
  
  // set parameters according to method in Ibe et al. (2013)
  void compute_spectra_MB_1loop();
  void compute_spectra_MB_2loop();
  
  // compute self energies using TSIL
  void compute_tsil();
  bool compute_tsil_iterative();
  
  // compute the iterative pole mass
  double iterative_ms_bar_mass(Data data, string particle);  
  
  double get_deltam();
  double get_deltam_2loop();
  
  double get_deltam2();
  double get_deltam2_2loop();
  
  double get_charged_mass();
  
  double get_double_charged_mass();
  
  double get_neutral_mass();
  
  double get_charged_mass_2loop();
  
  double get_neutral_mass_2loop();
  
  
};

#endif
