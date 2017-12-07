/*
 Mass Builder
 
 -- figures.cpp --
 
 Contains functions for generating the data for plots
 
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "figures_2.hpp"

using namespace std;

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

template <class T>
double get_deltam(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return spec.get_deltam() + spec.get_deltam_2loop();
}


template <class T>
double get_deltam2(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return spec.get_deltam2() + spec.get_deltam2_2loop();
}

template <class T>
double get_deltam_1loop(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	data.do_tsil_all = false;
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
  data.do_tsil_all = true;
	
  return spec.get_deltam() ;
}

template <class T>
double get_minus_deltam(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return -(spec.get_deltam() + spec.get_deltam_2loop());
}


template <class T>
double get_minus_deltam2(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return -(spec.get_deltam2() + spec.get_deltam2_2loop());
}


template <class S>
double find_min_Q(gsl_function F, double lower, double upper)
{
				
	int status;
	int iteration = 0, max_iteration = 150;
	
	cout << "lower , upper = " << lower << " , " << upper << endl;
	
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	
	double min = 0.5*(upper-lower)+lower;
	
	gsl_min_fminimizer_set (s, &F, min , lower, upper);
	
	do
	{
		iteration++;
		status = gsl_min_fminimizer_iterate (s);
		min = gsl_min_fminimizer_x_minimum (s);
		lower = gsl_min_fminimizer_x_lower (s);
		upper = gsl_min_fminimizer_x_upper (s);
		cout << "lower = " << lower << " upper = " << upper << endl;
		
		status = gsl_min_test_interval (lower, upper, 0.00001, 0.00001);
		
	}
	while (status == -2 && iteration < max_iteration);
	
	if (iteration == max_iteration)
	{
		cout << "could not determine minimum/maximum mass splitting" << endl;
	}
	
	gsl_min_fminimizer_free (s);
	
	return min;
	
}



template <class T>
Limits find_uncertainties2(Data data)
{
	
	Limits limits;
	gsl_function F;
	F.params = &data;
	
	bool do_search = false;
	
	double upper = 2.0*data.mt;
	double lower = 0.5*data.mt;
	
	// if within this specific range then the maximum and minimum
	// can be found within the Q range, not just at the end points
	
	if ( (data.MChi < 330) && (data.MChi > 200) && do_search )
	{
	  F.function = &get_deltam2<T>;
	  
	  limits.lower_Q = find_min_Q<T>(F, lower, upper);
	  limits.lower = get_deltam2<T>(limits.lower_Q,&data);
	  
		// find the maximum mass splitting
	  F.function = &get_minus_deltam2<T>;
	  
	  limits.upper_Q = find_min_Q<T>(F, lower, upper);
	  limits.upper = get_deltam2<T>(limits.upper_Q,&data);
	}
	else
	{
		double deltam_1 = get_deltam2<T>(upper,&data);
		double deltam_2 = get_deltam2<T>(lower,&data);
		if (deltam_1 < deltam_2)
		{
			limits.lower_Q = upper;
			limits.upper_Q = lower;
			limits.lower = deltam_1;
			limits.upper = deltam_2;
		}
		else
		{
			limits.lower_Q = lower;
			limits.upper_Q = upper;
			limits.lower = deltam_2;
			limits.upper = deltam_1;
		}
	}
	
	return limits;
	
}

template <class T>
Limits find_uncertainties(Data data)
{
	Limits limits;
	gsl_function F;
	F.params = &data;
	
	double upper = 2.0*data.mt;
	double lower = 0.5*data.mt;
	
	bool do_search = false; 
	 
	// if within this specific range then the maximum and minimum
	// can be found within the Q range, not just at the end points
	
	if ( (data.MChi < 300) && (data.MChi > 200) && do_search )
	{
	  F.function = &get_deltam<T>;
	  
	  limits.lower_Q = find_min_Q<T>(F, lower, upper);
	  limits.lower = get_deltam<T>(limits.lower_Q,&data);
	  
		// find the maximum mass splitting
	  F.function = &get_minus_deltam<T>;
	  
	  limits.upper_Q = find_min_Q<T>(F, lower, upper);
	  limits.upper = get_deltam<T>(limits.upper_Q,&data);
	}
	else
	{
		double deltam_1 = get_deltam<T>(upper,&data);
		double deltam_2 = get_deltam<T>(lower,&data);
		if (deltam_1 < deltam_2)
		{
			limits.lower_Q = upper;
			limits.upper_Q = lower;
			limits.lower = deltam_1;
			limits.upper = deltam_2;
		}
		else
		{
			limits.lower_Q = lower;
			limits.upper_Q = upper;
			limits.lower = deltam_2;
			limits.upper = deltam_1;
		}
	}
	
	return limits;
	
}



template <class T>
void Figures_2<T>::two_loop_plots(Data data, string group)
{ 
	int j = 1;
  if (group=="MDM" || group=="MDM_AB")
  {
		j = 2;
	}
	// output the maximum and minimum mass splittings (1 and 2-loop)
	string c_uncertainties = "mass_splittings/data2/uncertainties_" + group + ".txt";
  const char *unc_file = c_uncertainties.c_str();		
  // output the mass splitting at Q = mt  		
	string c_file_5 = "mass_splittings/data2/deltam_1loop_" + group + ".txt";
  const char *file_5 = c_file_5.c_str();	
	string c_file_6 = "mass_splittings/data2/deltam_2loop_" + group + ".txt";
  const char *file_6 = c_file_6.c_str();	 
  string c_file_9 = "mass_splittings/data2/deltam2_1loop_" + group + ".txt";
  const char *file_9 = c_file_9.c_str();	   
	string c_file_10 = "mass_splittings/data2/deltam2_2loop_" + group + ".txt";
  const char *file_10 = c_file_10.c_str();	
	// output the upper and lower values for the lifetimes
  string c_decay = "mass_splittings/data2/decays_" + group + ".txt";
  const char *decay_file = c_decay.c_str();
  // output the BRs and lifetimes for Q = mt at 1-loop
  string c_decay_BR = "mass_splittings/data2/decays_BR_1loop_" + group + ".txt";
  const char *decay_BR_file = c_decay_BR.c_str();
  // output the BRs and lifetimes for Q = mt at 2-loop
  string c_decay_BR_2loop = "mass_splittings/data2/decays_BR_2loop_" + group + ".txt";
  const char *decay_BR_file_2loop = c_decay_BR_2loop.c_str();
  
  string c_decay2_BR = "mass_splittings/data2/decays2_BR_1loop_" + group + ".txt";
  const char *decay2_BR_file = c_decay2_BR.c_str();
  // output the BRs and lifetimes for Q = mt at 2-loop
  
  string c_decay2_BR_2loop = "mass_splittings/data2/decays2_BR_2loop_" + group + ".txt";
  const char *decay2_BR_file_2loop = c_decay2_BR_2loop.c_str();
  
	
	ofstream unc_out;
	unc_out.open (unc_file);
  ofstream decays_out;
  decays_out.open (decay_file);  
  ofstream decays_BR_1loop;
  decays_BR_1loop.open (decay_BR_file);  
  ofstream decays_BR_2loop;
  decays_BR_2loop.open (decay_BR_file_2loop);  
  
  ofstream decays2_BR_1loop;
  ofstream decays2_BR_2loop;
  if (group=="MDM" || group=="MDM_AB")
  {
	  decays2_BR_1loop.open (decay2_BR_file);  
	  decays2_BR_2loop.open (decay2_BR_file_2loop);  
	}
  
  ofstream deltam_1loop;
  deltam_1loop.open (file_5);
  ofstream deltam_2loop;
  deltam_2loop.open (file_6);  
  
  // if MDM is the model then open these files
  ofstream deltam2_1loop;
  ofstream deltam2_2loop;
  if (group=="MDM" || group == "MDM_AB")
  {
		deltam2_1loop.open(file_9);
		deltam2_2loop.open(file_10);
	}
  
  // done opening files
  
	double mt = data.mt;
	double upper = 2.0*mt;
	double lower = 0.5*mt;
	
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(100.0L);
  //long double logMin = log10(1.0L);
  
  std::vector<double> Q(2);
  
  // number of points to plot
  int pts = 15;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
   
		Decays decays(data,j);
   
		unc_out << M ;
	  decays_out << M ;
	  decays_BR_1loop << M ;
	  decays_BR_2loop << M ;
	  
	  
	  
	  // get 1-loop uncertainties and decays
	  
	  // get decays and mass splittings at Q = mt
	  data.Q = mt;
	  
	  data.do_tsil_all = false;
	  T spec_1loop(data);
	  spec_1loop.compute_spectra_flexiblesusy();
	  spec_1loop.compute_tsil();
	    
	  deltam_1loop << spec_1loop.get_neutral_mass() ;
	  deltam_1loop << " " << spec_1loop.get_deltam();
	  
	  decays.calc_lifetime(spec_1loop.get_deltam(),decays_BR_1loop,0);
	    			
		if (group=="MDM" || group=="MDM_AB")
		{
			
			decays2_BR_1loop << M ;
			deltam2_1loop << spec_1loop.get_neutral_mass() ;
			deltam2_1loop << " " << spec_1loop.get_deltam2();
			decays.calc_lifetime(spec_1loop.get_deltam2(),decays2_BR_1loop,1);
		}
	  
	  // get uncertainties at 1-loop
	  data.Q = upper;
	  T spec_1loop_upper(data);
	  spec_1loop_upper.compute_spectra_flexiblesusy();
	  spec_1loop_upper.compute_tsil();
	  double deltam_1 = spec_1loop_upper.get_deltam();
	  data.Q = lower;
	  T spec_1loop_lower(data);
	  spec_1loop_lower.compute_spectra_flexiblesusy();
	  spec_1loop_lower.compute_tsil();
	  double deltam_2 = spec_1loop_lower.get_deltam();
	  
		if (deltam_1 < deltam_2)
		{
			unc_out << " " << deltam_1;
			unc_out << " " << deltam_2;
			decays_out << " " << decays.calc_lifetime(deltam_1,0);
			decays_out << " " << decays.calc_lifetime(deltam_2,0);
		}
		else
		{
			unc_out << " " << deltam_2;
			unc_out << " " << deltam_1;
			decays_out << " " << decays.calc_lifetime(deltam_2,0);
			decays_out << " " << decays.calc_lifetime(deltam_1,0);
		}
	  
	  if (group=="MDM" || group=="MDM_AB")
	  {
			deltam_1 = spec_1loop_upper.get_deltam2();
			deltam_2 = spec_1loop_lower.get_deltam2();
			if (deltam_1 < deltam_2)
			{
				unc_out << " " << deltam_1;
				unc_out << " " << deltam_2;
				decays_out << " " << decays.calc_lifetime(deltam_1,1);
				decays_out << " " << decays.calc_lifetime(deltam_2,1);
			}
			else
			{
				unc_out << " " << deltam_2;
				unc_out << " " << deltam_1;
				decays_out << " " << decays.calc_lifetime(deltam_2,1);
				decays_out << " " << decays.calc_lifetime(deltam_1,1);
			}
		}
	  
	  
	  
	  // get two-loop uncertainties
	  data.do_tsil_all = true;
	  
	  Limits limits = find_uncertainties<T>(data);
	  
	  unc_out << " " << limits.lower << " " << limits.upper;
	  decays_out << " " << decays.calc_lifetime(limits.lower,0);
	  decays_out << " " << decays.calc_lifetime(limits.upper,0);
	  
	  
	  if (group=="MDM" || group == "MDM_AB")
	  {
			
			Limits limits2 = find_uncertainties2<T>(data);
			
			unc_out << " " << limits2.lower << " " << limits2.upper;
			decays_out << " " << decays.calc_lifetime(limits2.lower,1);
			decays_out << " " << decays.calc_lifetime(limits2.upper,1);
			
		}
	    
		
		// get two-loop mass splittings and decays for Q = mt
		
		data.Q = mt;
		T spec_2loop(data);
	  spec_2loop.compute_spectra_flexiblesusy();
	  spec_2loop.compute_tsil();
	  
	  deltam_2loop << spec_2loop.get_neutral_mass() ;
	  deltam_2loop << " " << spec_2loop.get_deltam() + spec_2loop.get_deltam_2loop();
		decays.calc_lifetime(spec_2loop.get_deltam() + spec_2loop.get_deltam_2loop() ,decays_BR_2loop,0);
		
		
		if (group=="MDM" || group == "MDM_AB")
	  {
			decays2_BR_2loop << M ;
			decays.calc_lifetime(spec_2loop.get_deltam2() + spec_2loop.get_deltam2_2loop() ,decays2_BR_2loop,1);
			
			deltam2_2loop << spec_2loop.get_neutral_mass() ;
			
			deltam2_2loop << " " <<  spec_2loop.get_deltam2() + spec_2loop.get_deltam2_2loop();
			
			deltam2_1loop << endl;
			deltam2_2loop << endl;		
			
		}
		
		deltam_1loop << endl;
		deltam_2loop << endl;
	  unc_out << endl;
	  decays_out << endl;
		
  }
 
	deltam_1loop.close();
	deltam_2loop.close();
	deltam2_1loop.close();
	deltam2_2loop.close();
  unc_out.close();
  decays_out.close();
  decays_BR_1loop.close();
  decays_BR_2loop.close();
  
}



template <class T>
void Figures_2<T>::two_loop_plots_uncertainties(Data data, string group)
{ 
	double Pi = PI;
	// output the maximum and minimum mass splittings (1 and 2-loop)
	string c_uncertainties = "mass_splittings/data2/estimated_uncertainties_" + group + ".txt";
  const char *unc_file = c_uncertainties.c_str();		
  // output the mass splitting at Q = mt  	
  
  string c_decay = "mass_splittings/data2/estimated_decays_" + group + ".txt";
  const char *decay_file = c_decay.c_str();
  
  ofstream decays_out;
  decays_out.open (decay_file); 
  
  ofstream unc_out;
  unc_out.open (unc_file); 
  
	
  // done opening files
  
	double mt = data.mt;
	
	double upper = 2.0*mt;
	double lower = 0.5*mt;
	
  long double logMax = log10(1.0e4L);
  long double logMin = log10(1.0L);
  
  //long double logMax = log10(300.0L);
  //long double logMin = log10(200.0L);
  
  int j = 1;
  if (group=="MDM" || group=="MDM_AB")
  {
		j = 2;
	}
		
  // number of points to plot
  int pts = 20;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
   
		Decays decays(data,j);
   
		unc_out << M ;
	  decays_out << M ;
	  
	  // get spectrum at upper Q
	  data.Q = mt;
	  data.do_tsil_all = true;
	  T spec_upper(data);
	  spec_upper.compute_spectra_flexiblesusy();
	    
	  //double cw =  spec_upper.data.mw/spec_upper.data.mz;
	  double cw =  data.mw/data.mz;
    double cw2 =  pow(cw,2);
    double sw =  pow(1.-cw2,0.5);
	  double e =  pow(4*Pi*spec_upper.data.alpha,0.5);
    double g2 =  e/sw;
    double alpha2 = pow(g2,2)/(4*Pi);
	  
	  double error_upper = pow(alpha2/(4*Pi) , 2) * Pi * mt; // using top pole mass

	  cout << "error = " << error_upper*1000 << " MeV" << endl;
	  
	  double error_lower = error_upper;
	  
	  data.do_tsil_all = false;
	  
		// get uncertainties at 1-loop
	  data.Q = upper;
	  T spec_1loop_upper(data);
	  spec_1loop_upper.compute_spectra_flexiblesusy();
	  spec_1loop_upper.compute_tsil();
	  double deltam_1 = spec_1loop_upper.get_deltam();
	  
	  data.Q = lower;
	  T spec_1loop_lower(data);
	  spec_1loop_lower.compute_spectra_flexiblesusy();
	  spec_1loop_lower.compute_tsil();
	  double deltam_2 = spec_1loop_lower.get_deltam();
	 // double error_2 = pow( pow(spec_1loop_lower.data.g2,2)/(16*Pi*Pi) , 2) * Pi * mt;
	  
	  
	  
	  
		if (deltam_1 < deltam_2)
		{
			unc_out << " " << deltam_1-error_upper;
			unc_out << " " << deltam_2+error_lower;
			decays_out << " " << decays.calc_lifetime(deltam_1-error_upper,0);
			decays_out << " " << decays.calc_lifetime(deltam_2+error_lower,0);
		}
		else
		{
			unc_out << " " << deltam_2-error_lower;
			unc_out << " " << deltam_1+error_upper;
			decays_out << " " << decays.calc_lifetime(deltam_2-error_upper,0);
			decays_out << " " << decays.calc_lifetime(deltam_1+error_lower,0);
		}
	  
	  if (group=="MDM" || group=="MDM_AB")
	  {
			deltam_1 = spec_1loop_upper.get_deltam2();
			deltam_2 = spec_1loop_lower.get_deltam2();
			if (deltam_1 < deltam_2)
			{
				unc_out << " " << deltam_1-4.0*error_upper;
				unc_out << " " << deltam_2+4.0*error_lower;
				decays_out << " " << decays.calc_lifetime(deltam_1-4.0*error_upper,1);
				decays_out << " " << decays.calc_lifetime(deltam_2+4.0*error_lower,1);
			}
			else
			{
				unc_out << " " << deltam_2-4.0*error_lower;
				unc_out << " " << deltam_1+4.0*error_upper;
				decays_out << " " << decays.calc_lifetime(deltam_2-4.0*error_lower,1);
				decays_out << " " << decays.calc_lifetime(deltam_1+4.0*error_upper,1);
			}
		}
		
	  unc_out << endl;
	  decays_out << endl;
		
  }
 
	decays_out.close();
	unc_out.close();
  
}



template <class T>
void Figures_2<T>::test_MDM(Data data)
{ 
	
  
	
	ofstream out_file;
	out_file.open ("mass_splittings/data2/MDM_test.txt");
  
  
  // done opening files
  
	double mt = data.mt;
	double upper = 2.0*mt;
	double lower = 0.5*mt;
	
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(100);
  
  
  
  // number of points to plot
  int pts = 20;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
   
		
		out_file << M ;
	  
	  
	  
	  // get 1-loop uncertainties and decays
	  
	  // get decays and mass splittings at Q = mt
	  data.Q = mt;
	  
	  data.do_tsil_all = false;
	  T spec_1loop(data);
	  spec_1loop.compute_spectra_flexiblesusy(1,true);
	  spec_1loop.compute_tsil();
	    
	  out_file << " " << spec_1loop.get_deltam();
	  
	  
	  data.do_tsil_all = true;
	  
		T spec_2loop(data);
	  spec_2loop.compute_spectra_flexiblesusy(1,true);
	  spec_2loop.compute_tsil();
	  
	  
	  out_file << " " << spec_2loop.get_deltam() + spec_2loop.get_deltam_2loop();
		
		out_file << endl;
		
  }
 
	out_file.close();
  
}



template class Figures_2<EW_triplet_spectrum>;
template class Figures_2<MDM_spectrum>;
