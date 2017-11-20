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


template <class T>
void get_deltam(Data data, string parameter, long double central, long double range)
{
	Decays decays(data);
	cout << central-range << " < " << parameter << " < " << central+range << "  -> ";
	
	data.set_parameter(parameter,central-range);
	
	T specl(data);
  specl.compute_spectra_flexiblesusy(1);
	specl.compute_tsil();
	
	double deltam_upper = specl.get_deltam() + specl.get_deltam_2loop();
	double lifetime_upper = decays.calc_lifetime(deltam_upper);
	
  cout << deltam_upper  ;
  cout << "  < delta M < ";  
  
	data.set_parameter(parameter,central+range);
	
	T specu(data);
  specu.compute_spectra_flexiblesusy(1);
	specu.compute_tsil();  
	
	double deltam_lower = specu.get_deltam() + specu.get_deltam_2loop();
	double lifetime_lower = decays.calc_lifetime(deltam_lower);
	
	cout <<  deltam_lower;
	
	cout << " , range = " << 1000*abs(deltam_upper - deltam_lower) << " MeV " << endl;
	cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower))<< "%" << endl;
	
}





template <class T>
void uncertainties(Data data)
{
// vary SM parameters to see range of deltaM obtained
	
	// determine range in alpha for 1/alpha 127.940(42)
	long double upper = 1.0L/(1.0L/data.alpha + .042L);
	long double alpha_range = abs(upper-data.alpha);
	
	// determine range and central value for Q
	long double Q_upper = 2.0L*data.mt;
	long double Q_lower = 0.5L*data.mt;
	long double Q_central = 0.5*(Q_upper+Q_lower);
	long double Q_range = Q_upper - Q_central;
	
	// test effect of light quark masses
	
	Decays decays(data);
	T spec_no_ql_mass(data);
  spec_no_ql_mass.compute_spectra_flexiblesusy(1,true);
	spec_no_ql_mass.compute_tsil();
	double deltam_no_ql_mass = spec_no_ql_mass.get_deltam() + spec_no_ql_mass.get_deltam_2loop();
	double lifetime_upper = decays.calc_lifetime(deltam_no_ql_mass);
	
	T spec_ql_mass(data);
  spec_ql_mass.compute_spectra_flexiblesusy(1);
	spec_ql_mass.compute_tsil();
	double deltam_ql_mass = spec_ql_mass.get_deltam() + spec_ql_mass.get_deltam_2loop();	
	double lifetime_lower = decays.calc_lifetime(deltam_ql_mass);

	cout << " quark masses on -> delta M = " << deltam_ql_mass;
	cout << " quark masses zero -> delta M = " << deltam_no_ql_mass;
	cout << " , range = ";
	cout << abs(deltam_ql_mass - deltam_no_ql_mass)*1000;
	cout << " MeV" << endl;
	
	cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower))<< "%" << endl;
	

  get_deltam<T>(data, "mw" , data.mw, 0.015 );
  get_deltam<T>(data, "mz" , data.mz, 0.0021);
	get_deltam<T>(data, "alpha" , data.alpha, alpha_range);
	get_deltam<T>(data, "Q" , Q_central , Q_range);
  get_deltam<T>(data, "mh" , data.mh, 1.6L);
  get_deltam<T>(data, "mt" , data.mt, 2.28L);
  get_deltam<T>(data, "md" , data.md, 0.96e-3);
	
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
