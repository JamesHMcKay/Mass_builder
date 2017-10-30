/*
 Mass Builder
 
 -- splittings.cpp --

 Master routine for calling plotting functions for mass
 splitting results
 
 */

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"
#include "figures.hpp"
#include "mssm_spectrum.hpp"
#include "ew_triplet_spectrum.hpp"
#include "mdm_spectrum.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Data data(options);
  
  
  
  MDM_spectrum spec(data);
  spec.compute_spectra_flexiblesusy();
  
  data.alpha = spec.data.alpha;
  
  MDM_spectrum spec2(data);
  spec2.compute_spectra_MB_2loop();
  //spec.compute_tsil();
	/*
  cout << "--- explicit --- " << endl;
  cout << "1-loop mass splitting = " <<  spec.get_deltam() << endl;
  cout << "2-loop mass splitting = " <<  spec.get_deltam_2loop() + spec.get_deltam() << endl;
  cout << "1-loop mass splitting 2 = " <<  spec.get_deltam2() << endl;
  cout << "2-loop mass splitting 2 = " <<  spec.get_deltam2_2loop() + spec.get_deltam2() << endl;
	
  EW_triplet_spectrum spec2(data);
  spec2.compute_spectra_flexiblesusy();
  spec2.compute_tsil_iterative();
  
  cout << "--- iterative --- " << endl;
  cout << "1-loop mass splitting = " <<  spec2.get_deltam() << endl;
  cout << "2-loop mass splitting = " <<  spec2.get_deltam_2loop() + spec2.get_deltam()<< endl;
  */
  
  //Figures<EW_triplet_spectrum> fig;
  Figures<MDM_spectrum> fig;
  
  //fig.plot_M(data);
  
  //fig.plot_2loop_uncertainties(data,false);
  
  //fig.plot_M_2loop_explicit(data);
  
  //fig.plot_M_2loop_iterative(data);
  
  //fig.plot_M(data);
  
  //fig.plot_decays(data);
  
  //fig.plot_M_flexiblesusy(data);
  //fig.plot_M_flexiblesusy_2loop(data,"MDM",false);
  
  //fig.plot_deltam_2loop(data);
  
  //fig.plot_uncertainties(data);
  
	//fig.test(data);
  
  return 0;
}
