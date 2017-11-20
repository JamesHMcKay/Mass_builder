/*
 Mass Builder
 
 -- figures.cpp --
 
 Contains functions for generating the data for plots
 
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "figures.hpp"

using namespace std;

template <class T>
bool check_positive_deltam(Data data, double scale)
{
	data.Q = scale;
	
	T spec(data);
	bool error = spec.compute_spectra_flexiblesusy();
	
	if (error) {return true;}
  spec.compute_tsil_iterative();
	
  double deltam =  spec.get_deltam();
  if (deltam < 0 )
  {
		return true;
	}
	return false;
}
	
// find upper limit on Q before iteration fails
template <class T>
double get_bad_Q(Data data, double lower, double upper)
{
	bool error_a, error_b, error_c;
	
	// search over Q in log space
	double max = log10(upper);
	double min = log10(lower);
	
	double a = min;
	double b = max;
	
	error_a = check_positive_deltam<T>(data,pow(10,a));

	error_b = check_positive_deltam<T>(data,pow(10,b));
	
	if (error_b && !error_a)
	{
		while (abs(a-b)>0.0000001)
		{
			double c = 0.5*(a-b) + b ;
			error_c = check_positive_deltam<T>(data,pow(10,c));
			
			if (error_c)
			{
				b = c;
			}
			else
			{
				a = c;
			}
		}
	}
	
	return b;	
}


template <class T>
double get_mc(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil_iterative();
	
  return spec.get_charged_mass();
}

template <class T>
double get_mn(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil_iterative();
	
  return spec.get_neutral_mass();
}

template <class T>
double get_mn_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return spec.get_neutral_mass();
}

template <class T>
double get_mc_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return spec.get_charged_mass();
}

template <class T>
double get_min_Q(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil_iterative();
	
  return spec.get_deltam();
}

template <class T>
double get_max_Q(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
  T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil_iterative();
	
  return -spec.get_deltam();
}

template <class T>
double get_max_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return -spec.get_deltam();
}

template <class T>
double get_min_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
  
	T spec(data);
	spec.compute_spectra_flexiblesusy();
  spec.compute_tsil();
	
  return spec.get_deltam();
}

template <class S>
double find_min_Q(gsl_function F, double lower, double upper,Data data, bool minimum)
{
				
	int status;
	int iteration = 0, max_iteration = 1000;
	
	if (minimum)
	{
		// carefully bracket the minimum (sometimes hard to find)
		
		long double logMax = log10(upper);
		long double logMin = log10(lower);
		
		// plot deltam over range of Q values
	
		ofstream myfile;
		myfile.open ("mass_splittings/data/deltam_Q.txt");	   
	    
		std::vector<double> deltam(99, 10), Qval(99, 10); 
	  bool error_encountered = false;
	  // number of points to plot
	  int pts = 100;
		
	  for (int i=1;i<pts+1;i++)
		{
	    double n = i*(logMax - logMin)/pts + logMin;		
			
			data.Q = pow(10.0L,n);
			
			double deltam_test = F.function(data.Q,&data);
			 
			S spec(data);
			bool error = spec.compute_spectra_flexiblesusy();
			
			if (error && !error_encountered) {error_encountered = true; break;};
			
			spec.compute_tsil_iterative();
			
			
			if (deltam_test < 0 )
			{
				deltam[i-1] = 100;
			}
			else
			{
				deltam[i-1] = deltam_test;
			}
			Qval[i-1] = data.Q;
			myfile << data.Q << " " << deltam_test << endl;		
		}
		myfile.close();
		std::vector<double>::iterator result = std::min_element(std::begin(deltam),std::end(deltam));
		
		unsigned int it = std::distance(std::begin(deltam), result);
	
		
		if (it == deltam.size()-1)
		{
			upper = Qval[it];
		}
		else
		{
			upper = Qval[it+1];
		}
		if (it == 0)
		{
			lower = Qval[it];
		}
		else
		{
			lower = Qval[it-1];
		}
			
			cout << "it = " << it << endl;
			cout << "minimum = " << deltam[it] << " found at " << Qval[it] << endl;
		
	}
	
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
void Figures<T>::plot_uncertainties(Data data)
{
	data.do_tsil_all = false;
	int xi_int = data.Xi;
	
	string c_uncertainties = "mass_splittings/data/uncertainties_" + to_string(xi_int) + ".txt";
  const char *unc_file = c_uncertainties.c_str();
  
  string c_decay_it = "mass_splittings/data/decays_iterative_" + to_string(xi_int) + ".txt";
  const char *decay_it_file = c_decay_it.c_str();
    
  string c_decay_ex = "mass_splittings/data/decays_explicit_" + to_string(xi_int) + ".txt";
  const char *decay_ex_file = c_decay_ex.c_str();
	
	ofstream myfile;
	myfile.open (unc_file);
	
	ofstream decays_iterative;
  decays_iterative.open (decay_it_file);
  
  ofstream decays_explicit;
  decays_explicit.open (decay_ex_file);  
  
  // set range of plot
  long double logMax = log10(10000);
  long double logMin = log10(10.0);
  
  double mt = 173.15;

	double improved_upper_limit;
	double improved_upper_limit_plus_1;  
	
	Decays decays(data);
  
  // number of points to plot
  int pts =  150;
  for (int i=1;i<pts+1;i++)
	{
		
    double n = i*(logMax - logMin)/pts + logMin;		
		
		data.MChi = pow(10.0L,n);
		data.P = pow(10.0L,n);
		
		decays_iterative << data.MChi;
		decays_explicit << data.MChi;
		
		double upper = 2.0*data.MChi;
		double lower = 0.5*mt;
		
		if ( data.MChi < mt )
		{
			lower = 0.5*data.MChi;
			
			upper = 2.0*mt;
		}
		
		gsl_function F;
		F.params = &data;
	
		// iterative calculation
	
		cout << "--- iterative --- " << endl;
	
		// find the minimum mass splitting
		
		
		if (data.MChi > 700 )
		{
			improved_upper_limit = pow( 10.0 , get_bad_Q<T>(data,lower, upper) - 0.1 ) ;
		  improved_upper_limit_plus_1 = pow( 10.0 , get_bad_Q<T>(data,lower, upper) ) ;
		}
		else
		{
			improved_upper_limit = pow( 10.0 , get_bad_Q<T>(data,lower, upper)) ;
		  improved_upper_limit_plus_1 = pow( 10.0 , get_bad_Q<T>(data,lower, upper) ) ;
		}
		
		cout << "data.MChi = " << data.MChi << endl;
		cout << "previous upper limt = " << upper << " improved upper limit = " << improved_upper_limit << endl;
		
	  F.function = &get_min_Q<T>;
	  
	  double minQ_iterative = find_min_Q<T>(F, lower, improved_upper_limit, data, true);
		
		cout << "minimum delta M " << get_min_Q<T>(minQ_iterative,&data) << " at Q = " << minQ_iterative << endl;
	
		decays_iterative << " " << decays.calc_lifetime(get_min_Q<T>(minQ_iterative,&data));
	
		// find the maximum mass splitting
	  F.function = &get_max_Q<T>;
		double maxQ_iterative = find_min_Q<T>(F, lower, improved_upper_limit_plus_1, data, false);
				
		cout << "maximum delta M " << -get_max_Q<T>(maxQ_iterative,&data) << " at Q = " << maxQ_iterative << endl;
	
		decays_iterative << " " << decays.calc_lifetime(-get_max_Q<T>(maxQ_iterative,&data));
	
		// non-iterative calculation
		
		cout << "--- non-iterative --- " << endl;
		
	  F.function = &get_min_Q_low<T>;
	  
	  double minQ_non_iterative = find_min_Q<T>(F, lower, upper, data, false);
		
		cout << "minimum delta M " << get_min_Q_low<T>(minQ_non_iterative,&data) << " at Q = " << minQ_non_iterative << endl;
	
		decays_explicit << " " << decays.calc_lifetime(get_min_Q_low<T>(minQ_non_iterative,&data));
	
		// find the maximum mass splitting
	  F.function = &get_max_Q_low<T>;
	  
	  double maxQ_non_iterative = find_min_Q<T>(F, lower, upper, data, false);
		cout << "maximum delta M " << -get_max_Q_low<T>(maxQ_non_iterative,&data) << " at Q = " << maxQ_non_iterative << endl;

		decays_explicit << " " << decays.calc_lifetime(-get_max_Q_low<T>(maxQ_non_iterative,&data));

		// plot with respect to tree-level mass
		myfile << data.MChi << " " << get_min_Q<T>(minQ_iterative,&data) << " " << -get_max_Q<T>(maxQ_iterative,&data);
		myfile << " " << get_min_Q_low<T>(minQ_non_iterative,&data) << " " << -get_max_Q_low<T>(maxQ_non_iterative,&data) << endl;
		
		
		decays_explicit << endl;
		decays_iterative << endl;
		
	}
	
	myfile.close();
	decays_explicit.close();
	decays_iterative.close();
	
}

template <class T>
void Figures<T>::plot_Q(Data data)
{
  ofstream myfile;
  myfile.open ("mass_splittings/data2/deltam_Q.txt");
  int pts = 5;
  int status = 0;
  double max_Q = 400; // (GeV)
  double min_Q = 50; // (GeV)
  for (int i = 0; i < pts+1 ; i++)
  {
    data.mt = 163.3;  
    data.Q = (i)*(max_Q - min_Q)/pts + min_Q;
    
    data.do_tsil_all = false;
    T mssm_1loop(data);
    //mssm_1loop.compute_spectra_flexiblesusy();
    mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    double delta_m_1 = mssm_1loop.get_deltam();
    
    data.do_tsil_all = true;
    T mssm_2loop(data);
    //mssm_2loop.compute_spectra_flexiblesusy();
    mssm_2loop.compute_spectra_MB_2loop();
    mssm_2loop.compute_tsil();
    
    double delta_m_2 = mssm_2loop.get_deltam_2loop() + mssm_2loop.get_deltam() ;
    
    myfile << data.Q << " " <<   delta_m_1 <<  " " << delta_m_2 << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  myfile.close();  
}

// used to replicate Figure in Ibe et al.
template <class T>
void Figures<T>::plot_M(Data data)
{
  ofstream myfile;
  myfile.open ("mass_splittings/data2/deltam_M_lightquarks.txt");
  int pts = 60;
  double n = 0;
  int status = 0;
  double max_M = 5000; // (GeV)
  double min_M = 90; // (GeV)
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    data.MChi = pow(10,n);
    data.P = data.MChi;
    data.do_tsil_all = false;
    T mssm_1loop(data);
    mssm_1loop.compute_spectra_flexiblesusy();
    
    //mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    double delta_m_1 = mssm_1loop.get_deltam();
    
    data.do_tsil_all = true;
    
    
    T mssm_2loop_1loopRGE(data);
    
    mssm_2loop_1loopRGE.compute_spectra_flexiblesusy(1);
    mssm_2loop_1loopRGE.compute_tsil();
    double delta_m_2_1 = mssm_2loop_1loopRGE.get_deltam_2loop()+ mssm_2loop_1loopRGE.get_deltam();
    
    T mssm_2loop_2loopRGE(data);
    
    mssm_2loop_2loopRGE.compute_spectra_flexiblesusy(2);
    mssm_2loop_2loopRGE.compute_tsil();
    double delta_m_2_2 = mssm_2loop_2loopRGE.get_deltam_2loop()+ mssm_2loop_2loopRGE.get_deltam();
    
    //myfile << data.MChi << " " << delta_m_1 <<  " " << delta_m_2;
    //myfile << " " << data.MChi << " " << delta_m_2FS << endl;
    
    myfile << mssm_2loop_1loopRGE.get_neutral_mass() << " " << delta_m_1 <<  " " << delta_m_2_1;
    myfile << " " << mssm_2loop_2loopRGE.get_neutral_mass() << " " << delta_m_2_2 << endl;
    
    status=(float(i)/pts)*100;
		cout<< "\r" << "computing mass splitting . . . " << status << "% complete ";
	  std::cout << std::flush;
  }
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  myfile.close();
}

template <class T>
void Figures<T>::plot_M_flexiblesusy(Data data)
{ 
	data.do_tsil_all = false;
	
	int xi_int = data.Xi;
	
	string c_deltam_it = "mass_splittings/data/mass_splittings_iterative_" + to_string(xi_int) + ".txt";
  const char *deltam_it_file = c_deltam_it.c_str();
  
	string c_deltam_ex = "mass_splittings/data/mass_splittings_explicit_" + to_string(xi_int) + ".txt";
  const char *deltam_ex_file = c_deltam_ex.c_str();
  
	string c_pole_mass_c_it = "mass_splittings/data/pole_mass_c_iterative_" + to_string(xi_int) + ".txt";
  const char *pole_mass_c_it_file = c_pole_mass_c_it.c_str();
  
	string c_pole_mass_n_it = "mass_splittings/data/pole_mass_n_iterative_" + to_string(xi_int) + ".txt";
  const char *pole_mass_n_it_file = c_pole_mass_n_it.c_str();
  
	string c_pole_mass_c_ex = "mass_splittings/data/pole_mass_c_explicit_" + to_string(xi_int) + ".txt";
  const char *pole_mass_c_ex_file = c_pole_mass_c_ex.c_str();
  
	string c_pole_mass_n_ex = "mass_splittings/data/pole_mass_n_explicit_" + to_string(xi_int) + ".txt";
  const char *pole_mass_n_ex_file = c_pole_mass_n_ex.c_str();          
	
	
  ofstream mass_splittings_iterative;
  mass_splittings_iterative.open (deltam_it_file);
  
  ofstream mass_splittings_explicit;
  mass_splittings_explicit.open (deltam_ex_file);  
 
  ofstream pole_mass_n_iterative;
  pole_mass_n_iterative.open (pole_mass_n_it_file);
 
	ofstream pole_mass_c_iterative;
  pole_mass_c_iterative.open (pole_mass_c_it_file);
  
  ofstream pole_mass_n_explicit;
  pole_mass_n_explicit.open (pole_mass_n_ex_file);
 
	ofstream pole_mass_c_explicit;
  pole_mass_c_explicit.open (pole_mass_c_ex_file);  
 
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 200;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
    
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    mass_splittings_iterative << M ;
    pole_mass_n_iterative << M ;
    pole_mass_c_iterative << M ;
    
    mass_splittings_explicit << M ;
    pole_mass_n_explicit << M ;
    pole_mass_c_explicit << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];
			T spec(data);
			spec.compute_spectra_flexiblesusy();
			
			T spec_iterative = spec;
      bool error = spec_iterative.compute_tsil_iterative();
      if (!error)
      {
							mass_splittings_iterative << " " << 0.0;			
							pole_mass_n_iterative << " " << 0.0;
							pole_mass_c_iterative << " " << 0.0;
			}
			else
			{
				mass_splittings_iterative << " " << spec_iterative.get_deltam();			
				pole_mass_n_iterative << " " << spec_iterative.get_neutral_mass();
				pole_mass_c_iterative << " " << spec_iterative.get_charged_mass();
			}
			T spec_explicit = spec;
      spec_explicit.compute_tsil();
			
			mass_splittings_explicit << " " << spec_explicit.get_deltam();
			
			pole_mass_n_explicit << " " << spec_explicit.get_neutral_mass();
			pole_mass_c_explicit << " " << spec_explicit.get_charged_mass();
			
		}
    
    mass_splittings_iterative << endl;
    pole_mass_n_iterative << endl;
    pole_mass_c_iterative << endl;
    
		mass_splittings_explicit << endl;
    pole_mass_n_explicit << endl;
    pole_mass_c_explicit << endl;    
  }
  mass_splittings_iterative.close();
  pole_mass_n_iterative.close();
  pole_mass_c_iterative.close();
  mass_splittings_explicit.close();
  pole_mass_n_explicit.close();
  pole_mass_c_explicit.close();  
  
}


template <class T>
void Figures<T>::plot_M_flexiblesusy_2loop(Data data, string group,bool do_iteration)
{ 
		
	string c_file_1 = "mass_splittings/data2/pole_mass_unc_2loop_n_" + group + ".txt";
  const char *file_1 = c_file_1.c_str();	
		
	string c_file_2 = "mass_splittings/data2/pole_mass_unc_2loop_c_" + group + ".txt";
  const char *file_2 = c_file_2.c_str();	
    		
	string c_file_3 = "mass_splittings/data2/pole_mass_unc_1loop_n_" + group + ".txt";
  const char *file_3 = c_file_3.c_str();	
		
	string c_file_4 = "mass_splittings/data2/pole_mass_unc_1loop_c_" + group + ".txt";
  const char *file_4 = c_file_4.c_str();	
    		
	string c_file_5 = "mass_splittings/data2/deltam_1loop_" + group + ".txt";
  const char *file_5 = c_file_5.c_str();	
  
	string c_file_6 = "mass_splittings/data2/deltam_2loop_" + group + ".txt";
  const char *file_6 = c_file_6.c_str();	 
  
	string c_file_7 = "mass_splittings/data2/decays_1loop_" + group + ".txt";
  const char *file_7 = c_file_7.c_str();	   
 
	string c_file_8 = "mass_splittings/data2/decays_2loop_" + group + ".txt";
  const char *file_8 = c_file_8.c_str();	     
  
  
  
  string c_file_9 = "mass_splittings/data2/deltam2_1loop_" + group + ".txt";
  const char *file_9 = c_file_9.c_str();	   
 
	string c_file_10 = "mass_splittings/data2/deltam2_2loop_" + group + ".txt";
  const char *file_10 = c_file_10.c_str();	
  
 
  
	
  ofstream pole_mass_n_2loop;
  pole_mass_n_2loop.open (file_1);
 
	ofstream pole_mass_c_2loop;
  pole_mass_c_2loop.open (file_2);
  
  ofstream pole_mass_n_2loop_it;
  ofstream pole_mass_c_2loop_it;
  if (do_iteration)
  {
		pole_mass_n_2loop_it.open ("mass_splittings/data2/pole_mass_unc_2loop_n_it.txt");
		pole_mass_c_2loop_it.open ("mass_splittings/data2/pole_mass_unc_2loop_c_it.txt");  
	}

	ofstream pole_mass_n_1loop;
  pole_mass_n_1loop.open (file_3);
 
	ofstream pole_mass_c_1loop;
  pole_mass_c_1loop.open (file_4);
  
  ofstream deltam_1loop;
  deltam_1loop.open (file_5);
  
  ofstream deltam_2loop;
  deltam_2loop.open (file_6);  
  
  ofstream deltam_2loop_it;
  if (do_iteration)
  {
		deltam_2loop_it.open ("mass_splittings/data2/deltam_2loop_it.txt");    
  }
  
  ofstream deltam2_1loop;
  ofstream deltam2_2loop;
  
  if (group=="MDM" || group == "MDM_AB")
  {
		deltam2_1loop.open(file_9);
		deltam2_2loop.open(file_10);
	}
	
	// output file for decay lifetimes (non-iterative only)
	
	ofstream decays_1loop;
	decays_1loop.open(file_7);

	ofstream decays_2loop;
	decays_2loop.open(file_8);	
 
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(100.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 10;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
    
    
    //Q[0] = 2.0 * M;
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 173.15;
    //Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    pole_mass_n_2loop << M ;
    pole_mass_c_2loop << M ;
    
    pole_mass_n_1loop << M ;
    pole_mass_c_1loop << M ;
    
    deltam_1loop << M ;
    deltam_2loop << M ;
    deltam2_1loop << M ;
    deltam2_2loop << M ;    
    decays_1loop << M ;
    decays_2loop << M ;
    
    if (do_iteration)
    {
			deltam_2loop_it << M ;
			pole_mass_n_2loop_it << M ;
			pole_mass_c_2loop_it << M ;
		}
    
    for (int i = 0; i < 4 ; i++)
    {
			data.Q = Q[i];
			data.do_tsil_all = true;
			T spec(data);
			
			spec.compute_spectra_flexiblesusy();
			
			
			// do explicit calculation
			data.do_tsil_all = false;
	  	T spec_1loop(data);
	  	spec_1loop.compute_spectra_flexiblesusy();
	    spec_1loop.compute_tsil();
	    deltam_1loop << " " << spec_1loop.get_deltam();
	    			
			
			if (group=="MDM" || group=="MDM_AB")
			{
				deltam2_1loop << " " << spec_1loop.get_deltam2();
			}
	    
	  	T spec_explicit = spec;
	    spec_explicit.compute_tsil();
	    Decays decays(data);
				
			pole_mass_n_2loop << " " << spec_explicit.get_neutral_mass_2loop();
			pole_mass_c_2loop << " " << spec_explicit.get_charged_mass_2loop();	
			pole_mass_n_1loop << " " << spec_explicit.get_neutral_mass();	
			pole_mass_c_1loop << " " << spec_explicit.get_charged_mass();
			
			deltam_2loop << " " << spec_explicit.get_deltam_2loop()+spec_explicit.get_deltam();
			
			decays_1loop << " " << decays.calc_lifetime(spec_explicit.get_deltam());
			decays_2loop << " " << decays.calc_lifetime(spec_explicit.get_deltam_2loop()+spec_explicit.get_deltam());
			
			if (group=="MDM" || group=="MDM_AB")
			{
				deltam2_2loop << " " << spec_explicit.get_deltam2_2loop()+spec_explicit.get_deltam2();
			}
			
			// do iterative calculation
			
			if (do_iteration) // && data.Q < 1000 && data.MChi > 1000)
			{
				T spec_iterative = spec;
	      bool error = spec_iterative.compute_tsil_iterative();
	      if (!error)
	      {
					deltam_2loop_it << " " << 0.0L;
					pole_mass_n_2loop_it << " " << 0.0L;
					pole_mass_c_2loop_it << " " << 0.0L;
				}
				else
				{
					deltam_2loop_it << " " << spec_iterative.get_deltam() + spec_iterative.get_deltam_2loop();
					pole_mass_n_2loop_it << " " << spec_iterative.get_neutral_mass_2loop();
					pole_mass_c_2loop_it << " " << spec_iterative.get_charged_mass_2loop();				
				}
			}
			
			
		}
    
    pole_mass_n_2loop << endl;
    pole_mass_c_2loop << endl;
    
    pole_mass_n_1loop << endl;    
		pole_mass_c_1loop << endl;
		
		deltam_1loop << endl;
		deltam_2loop << endl;
		deltam2_1loop << endl;
		deltam2_2loop << endl;		
		
		decays_1loop << endl;
		decays_2loop << endl;
		
		if (do_iteration)
		{
			deltam_2loop_it << endl;
			pole_mass_n_2loop_it << endl;
			pole_mass_c_2loop_it << endl;		
		}
		
  }
  
  pole_mass_n_2loop.close();
  pole_mass_c_2loop.close();
  pole_mass_n_2loop_it.close();
  pole_mass_c_2loop_it.close();  
  pole_mass_n_1loop.close();
  pole_mass_c_1loop.close();
  
  deltam_1loop.close();
  deltam_2loop.close();
  deltam2_1loop.close();
  deltam2_2loop.close();  
  deltam_2loop_it.close();
  decays_1loop.close();
  decays_2loop.close();
  
}


template <class T>
void Figures<T>::plot_deltam_2loop(Data data)
{ 
  ofstream deltam_2loop_it;
  deltam_2loop_it.open ("mass_splittings/data2/deltam_2loop_it.txt");    
 
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  bool error = false;
  // number of points to plot
  int pts = 30; // 50;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
    
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    
		deltam_2loop_it << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];
			T spec(data);
			spec.compute_spectra_flexiblesusy();
			
			T spec_iterative = spec;

			error = spec_iterative.compute_tsil_iterative();
			
      if (!error)
      {
				deltam_2loop_it << " " << 0.0L;
			}
			else
			{
				if (spec_iterative.get_deltam_2loop() < 1.0)
				{
					deltam_2loop_it << " " << spec_iterative.get_deltam() + spec_iterative.get_deltam_2loop();
				}
				else
				{
					deltam_2loop_it << " " << 0.0L;
				}
			}
			
		}
		
		deltam_2loop_it << endl;
		
  }
  
  deltam_2loop_it.close();
}

template class Figures<EW_triplet_spectrum>;
template class Figures<MDM_spectrum>;
