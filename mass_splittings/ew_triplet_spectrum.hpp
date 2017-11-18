#ifndef EW_TRIPLET_H
#define EW_TRIPLET_H

#include "data.hpp"


class EW_triplet_spectrum
{
private:
		
public:
  
  // Mass Builder data structure
  Data data;
  double MChi_pole;
   
  EW_triplet_spectrum(Data data) : data(data) {};
  
  EW_triplet_spectrum(){};
  
  // run FlexibleSUSY and update data structure with couplings and self energies
  bool compute_spectra_flexiblesusy(int loop_order = 1);
  
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
  
  double get_charged_mass();
  double get_charged_mass_2loop();
  
  double get_neutral_mass();
  double get_neutral_mass_2loop();
  
  double get_F11_der();
  double get_F12_der();
  
  double get_F11_der_it();
  double get_F12_der_it();
  
  // dummy functions to satisfy templated functions
  double get_deltam2(){return 0;};
  double get_deltam2_2loop(){return 0;};
    
};

#endif
