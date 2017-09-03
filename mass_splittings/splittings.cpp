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
  
  MSSM_spectrum spec(data);
  spec.compute_spectra_flexiblesusy();
  
  Figures<EW_triplet_spectrum> fig;
  
  fig.plot_M(data);
  
  //fig.plot_decays(data);
  
  //fig.plot_M_flexiblesusy_2loop(data);
  
  //fig.plot_M_flexiblesusy(data);
  
  //fig.plot_uncertainties(data);
  
	//fig.test(data);
  
  return 0;
}
