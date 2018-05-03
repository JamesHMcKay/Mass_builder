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

bool MDM = true;
int j = 0;

template <class T>
void get_deltam(Data data, string parameter, long double central, long double range)
{
	Decays decays(data,j);
	cout << central-range << " < " << parameter << " < " << central+range << "  -> ";
	
	data.set_parameter(parameter,central-range);
	
	T specl(data);
  specl.compute_spectra_flexiblesusy(1);
	specl.compute_tsil();
	
	double deltam_upper = specl.get_deltam() + specl.get_deltam_2loop();
	double lifetime_upper = decays.calc_lifetime(deltam_upper,0);
	
  cout << deltam_upper  ;
  cout << "  < delta M < ";  
  
	data.set_parameter(parameter,central+range);
	
	T specu(data);
  specu.compute_spectra_flexiblesusy(1);
	specu.compute_tsil();  
	
	double deltam_lower = specu.get_deltam() + specu.get_deltam_2loop();
	double lifetime_lower = decays.calc_lifetime(deltam_lower,0);
	
	if (MDM)
	{
		double deltam_upper2 = specl.get_deltam2() + specl.get_deltam2_2loop();
		double lifetime_upper2 = decays.calc_lifetime(deltam_upper2,1);	
		double deltam_lower2 = specu.get_deltam2() + specu.get_deltam2_2loop();
		double lifetime_lower2 = decays.calc_lifetime(deltam_lower2,1);	
		
		cout <<  deltam_lower;	
		cout << " , range = " << 1000*abs(deltam_upper - deltam_lower);
		cout << " ( " << 1000*abs(deltam_upper2 - deltam_lower2) << " ) " << " MeV " << endl;
		cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower));
		cout << " ( " << 100*abs(lifetime_upper2-lifetime_lower2)/ ( 0.5*(lifetime_upper2+lifetime_lower2)) << " ) ";
		cout << "%" << endl;
	}
	else
	{
		cout <<  deltam_lower;	
		cout << " , range = " << 1000*abs(deltam_upper - deltam_lower) << " MeV " << endl;
		cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower))<< "%" << endl;
	}
}





template <class T>
void uncertainties(Data data)
{
	
	cout.precision(3);
	// force use of top mass for renormalisation scale
	data.Q = data.mt;
	
	get_deltam<T>(data,  "ml"     ,  data.ml     , 0.12e-3     );
	get_deltam<T>(data,  "mb"     ,  data.mb     , 0.09        );
	
	
	
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
	
	Decays decays(data,j);
	T spec_no_ql_mass(data);
  spec_no_ql_mass.compute_spectra_flexiblesusy(1,true);
	spec_no_ql_mass.compute_tsil();
	double deltam_no_ql_mass = spec_no_ql_mass.get_deltam() + spec_no_ql_mass.get_deltam_2loop();
	double lifetime_upper = decays.calc_lifetime(deltam_no_ql_mass,0);
	
	T spec_ql_mass(data);
  spec_ql_mass.compute_spectra_flexiblesusy(1);
	spec_ql_mass.compute_tsil();
	double deltam_ql_mass = spec_ql_mass.get_deltam() + spec_ql_mass.get_deltam_2loop();	
	double lifetime_lower = decays.calc_lifetime(deltam_ql_mass,0);
	
	
	double lifetime_lower2 = 0;
	double deltam_ql_mass2 = 0;
	
	if (MDM)
	{
		double deltam_no_ql_mass2 = spec_no_ql_mass.get_deltam2() + spec_no_ql_mass.get_deltam2_2loop();
		double lifetime_upper2 = decays.calc_lifetime(deltam_no_ql_mass2,1);
		deltam_ql_mass2 = spec_ql_mass.get_deltam2() + spec_ql_mass.get_deltam2_2loop();	
		lifetime_lower2 = decays.calc_lifetime(deltam_ql_mass2,1);
		cout << " quark masses on -> delta M (2-loop) = " << deltam_ql_mass;
		cout << " ( " << deltam_ql_mass2 << " ), ";
		cout << " quark masses zero -> delta M = " << deltam_no_ql_mass;
		cout << " ( " << deltam_no_ql_mass2 << " ) " ;
		cout << " , range = ";
		cout << abs(deltam_ql_mass - deltam_no_ql_mass)*1000;
		cout << " ( " << abs(deltam_ql_mass2 - deltam_no_ql_mass2)*1000 << " ) " ;
		cout << " MeV" << endl;
		
		cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower));
		cout << " ( " << 100*abs(lifetime_upper2-lifetime_lower2)/ ( 0.5*(lifetime_upper2+lifetime_lower2)) << " ) ";
		cout << "%" << endl;
	}
	else
	{
		cout << " quark masses on -> delta M (2-loop) = " << deltam_ql_mass;
		cout << " quark masses zero -> delta M = " << deltam_no_ql_mass;
		cout << " , range = ";
		cout << abs(deltam_ql_mass - deltam_no_ql_mass)*1000;
		cout << " MeV" << endl;
		
		cout << " effect on lifetime " << 100*abs(lifetime_upper-lifetime_lower)/ ( 0.5*(lifetime_upper+lifetime_lower))<< "%" << endl;
	}

// one-loop mass splitting

	T spec_1loop(data);
	spec_1loop.data.do_tsil_all = false;
  spec_1loop.compute_spectra_flexiblesusy(1);
	spec_1loop.compute_tsil();
	double deltam_1loop = spec_1loop.get_deltam();
	double lifetime_1loop = decays.calc_lifetime(deltam_1loop,0);
	
	if (MDM)
	{
		double deltam_1loop2 = spec_1loop.get_deltam2();
		double lifetime_1loop2 = decays.calc_lifetime(deltam_1loop2,1);
		
		cout << "1-loop mass splitting = " << deltam_1loop;
		cout << " ( " << deltam_1loop2 << " ) ";
		cout << " effect on lifetime compared to two-loop ";
		cout << 100*abs(lifetime_1loop-lifetime_lower)/ ( 0.5*(lifetime_1loop+lifetime_lower));
		cout << " ( " << 100*abs(lifetime_1loop2-lifetime_lower2)/ ( 0.5*(lifetime_1loop2+lifetime_lower2)) << " ) ";
		cout << "%" << endl;
	 
	 }


	get_deltam<T>(data,  "alpha"  ,  data.alpha  , alpha_range );
	get_deltam<T>(data,  "mt"     ,  data.mt     , 2.28L       );
  get_deltam<T>(data,  "mh"     ,  data.mh     , 1.6L				 );
  get_deltam<T>(data,  "mw"     ,  data.mw     , 0.015 			 );
  get_deltam<T>(data,  "mz"     ,  data.mz     , 0.0021			 );
	get_deltam<T>(data,  "me"     ,  data.me     , 0.31e-11		 );
	get_deltam<T>(data,  "mm"     ,  data.mm     , 0.24e-10    );
	get_deltam<T>(data,  "ml"     ,  data.ml     , 0.12e-3     );
	get_deltam<T>(data,  "md"     ,  data.md     , 0.96e-3     );
	get_deltam<T>(data,  "mu"     ,  data.mu     , 0.46e-3     );
	get_deltam<T>(data,  "ms"     ,  data.ms     , 15e-3       );
	get_deltam<T>(data,  "mc"     ,  data.mc     , 0.075       );
	get_deltam<T>(data,  "mb"     ,  data.mb     , 0.09        );
	get_deltam<T>(data,  "alphaS" ,  data.alphaS , 0.0011      );
	get_deltam<T>(data,  "Q"      ,  Q_central   , Q_range     );
	
}


int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Data data(options);
  
  // produce values for the uncertainties table
  // wino model
  //j = 1; uncertainties<EW_triplet_spectrum>(data);
  // MDM model
  //j = 2; uncertainties<MDM_spectrum>(data);
  
  
  // Figures for first paper
  //Figures<EW_triplet_spectrum> fig;
  
  //fig.plot_M(data);
  
  //fig.plot_M_2loop_explicit(data);
  
  //fig.plot_M_2loop_iterative(data);
  
  //fig.plot_M(data);
  
  //fig.plot_decays(data);
  
  //fig.plot_M_flexiblesusy(data);
  //fig.plot_M_flexiblesusy_2loop(data,"MSSM",false);
  
  //fig.plot_deltam_2loop(data);
  
  //fig.plot_uncertainties(data);
  
  
  
	// Figures for second paper
	
	/*
	Figures_2<EW_triplet_spectrum> fig_2;
	fig_2.two_loop_plots(data, "MSSM");
	fig_2.two_loop_plots_uncertainties(data, "MSSM");
	*/
	

	/*
	Figures_2<MDM_spectrum> fig_2;	
	fig_2.two_loop_plots(data, "MDM");
	fig_2.two_loop_plots_uncertainties(data, "MDM");
	*/
	
  return 0;
}
