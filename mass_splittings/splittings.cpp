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
#include "figures_2.hpp"
#include "ew_triplet_spectrum.hpp"
#include "mdm_spectrum.hpp"

using namespace std;

/*
template <class T>
double get_deltam(Data data, )
{
	
	T spec(data);
  spec.compute_spectra_flexiblesusy(1);
  
	spec.compute_tsil();
	
  cout << "2-loop mass splitting = " <<  spec.get_deltam() + spec.get_deltam_2loop() <<  endl;  
	
	
}*/





template <class T>
void uncertainties(Data data)
{
// vary SM parameters to see range of deltaM obtained


  T spec(data);
  spec.compute_spectra_flexiblesusy(1);
  
	spec.compute_tsil();
	
	cout.precision(17);
	
  cout << "1-loop mass splitting = " <<  spec.get_deltam() <<  endl;  
  cout << "2-loop mass splitting = " <<  spec.get_deltam() + spec.get_deltam_2loop() <<  endl;  
  
  // want to print out upper, lower and range when varying parameter x
  
  
	

}






int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Data data(options);
  
  
  uncertainties<EW_triplet_spectrum>(data);
  
  //EW_triplet_spectrum spec2(data);
  //spec2.compute_spectra_flexiblesusy();
  
	//spec2.compute_tsil();
	
  //cout << "MDM mass splitting = " <<  spec2.get_deltam() << " " << spec2.get_deltam_2loop() <<  endl;  
  
  /*

  spec.compute_tsil();
	
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
  
  Figures<EW_triplet_spectrum> fig;
  //Figures<MDM_spectrum> fig;
  
  //fig.plot_M(data);
  
  //fig.plot_M_2loop_explicit(data);
  
  //fig.plot_M_2loop_iterative(data);
  
  //fig.plot_M(data);
  
  //fig.plot_decays(data);
  
  //fig.plot_M_flexiblesusy(data);
  //fig.plot_M_flexiblesusy_2loop(data,"MSSM",false);
  
  //fig.plot_deltam_2loop(data);
  
  //fig.plot_uncertainties(data);
  
	//fig.test(data);
	
	Figures_2<EW_triplet_spectrum> fig_2;
	//Figures_2<MDM_spectrum> fig_2;
	
	//fig_2.two_loop_plots(data, "MDM");  
	
  return 0;
}
