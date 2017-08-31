/*
 Mass Builder
 
 -- figures.cpp --
 
 Contains functions for generating the data for plots
 
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "figures.hpp"

using namespace std;

bool check_positive_deltam(Data data, double scale)
{
	data.Q = scale;
	
	MSSM mssm(data);
	bool error = mssm.compute_spectra_flexiblesusy();
	
	if (error) {return true;}
	//mssm.compute_spectra_MB_1loop();
  mssm.compute_tsil_iterative();
	
	
	
  double deltam =  (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
  if (deltam < 0 )
  {
		return true;
	}
	return false;
}
	
// find upper limit on Q before iteration fails
double get_bad_Q(Data data, double lower, double upper)
{
	bool error_a, error_b, error_c;
	
	// search over Q in log space
	double max = log10(upper);
	double min = log10(lower);
	
	double a = min;
	double b = max;
	
	error_a = check_positive_deltam(data,pow(10,a));

	error_b = check_positive_deltam(data,pow(10,b));
	
	if (error_b && !error_a)
	{
		while (abs(a-b)>0.0000001)
		{
			double c = 0.5*(a-b) + b ;
			error_c = check_positive_deltam(data,pow(10,c));
			
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



double get_min_Q(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  mssm.compute_tsil_iterative();
	
  return (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double get_max_Q(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
  MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  mssm.compute_tsil_iterative();
	
  return -(mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}


double get_max_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  mssm.compute_tsil();
	
  return -(mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double get_min_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
  
	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  mssm.compute_tsil();
	
  return (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double find_min_Q(gsl_function F, double lower, double upper,Data data, bool minimum)
{
				
	int status;
	int iteration = 0, max_iteration = 1000;
	
	//Data data;
	//data.MChi = M; // (GeV)
	
	
	//Data data = *(Data*)F.params;
	
	cout << "M = " << data.MChi << endl;
	
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
			 
			MSSM mssm(data);
			bool error = mssm.compute_spectra_flexiblesusy();
			//cout << "error = " << error << " Q = " << data.Q << endl;
			
			if (error && !error_encountered) {error_encountered = true; break;};
			
			mssm.compute_tsil_iterative();
			
			
			if (/*mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]*/ deltam_test < 0 )
			{
				deltam[i-1] = 100;
			}
			else
			{
				deltam[i-1] = deltam_test;//mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"];
			}
			Qval[i-1] = data.Q;
			myfile << data.Q << " " << /*mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]*/ deltam_test << endl;		
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



void Figures::test(Data data)
{
		
  double mt = 173.15;

	double improved_upper_limit;

	double upper = 2.0*data.MChi;
	double lower = 0.5*mt;
	
	if ( data.MChi < mt )
	{
		lower = 0.5*data.MChi;
		
		upper = 2.0*mt;
	}
	
	gsl_function F;
	F.params = &data;
	// find the minimum mass splitting
	
	if (data.MChi > 700 )
	{
		improved_upper_limit = pow( 10.0 , get_bad_Q(data,lower, upper)-0.01  ) ;
	}
	else
	{
		improved_upper_limit = pow( 10.0 , get_bad_Q(data,lower, upper)) ;
	}
	
	cout << "data.MChi = " << data.MChi << endl;
	cout << "previous upper limt = " << upper << " improved upper limit = " << improved_upper_limit << endl;
	
  F.function = &get_min_Q_low;
  
  double minQ_non_iterative = find_min_Q(F, lower, improved_upper_limit, data,true);
	
	cout << "minimum delta M " << get_min_Q_low(minQ_non_iterative,&data) << " at Q = " << minQ_non_iterative << endl;

	data.Q =  minQ_non_iterative ;

	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	
	MSSM mssm_explicit = mssm;
	mssm_explicit.compute_tsil();
	
	cout<< "actual mass splitting at " <<  minQ_non_iterative << " = "  << " " << mssm_explicit.data.SE_1["F12_g1"] - mssm_explicit.data.SE_1["F11_g1"] << endl;

}


void Figures::plot_uncertainties(Data data)
{
	data.do_tsil_all = false;
	ofstream myfile;
	myfile.open ("mass_splittings/data/uncertainties.txt");
	
	ofstream decays_iterative;
  decays_iterative.open ("mass_splittings/data/decays_iterative.txt");
  
  ofstream decays_explicit;
  decays_explicit.open ("mass_splittings/data/decays_explicit.txt");  
  
  // set range of plot
  long double logMax = log10(10000);
  long double logMin = log10(10.0);
  
  double mt = 173.15;

	double improved_upper_limit;
	double improved_upper_limit_plus_1;  
	
	Decays decays(data);
  
  // number of points to plot
  int pts = 150;
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
			improved_upper_limit = pow( 10.0 , get_bad_Q(data,lower, upper) - 0.1 ) ;
		  improved_upper_limit_plus_1 = pow( 10.0 , get_bad_Q(data,lower, upper) ) ;
		}
		else
		{
			improved_upper_limit = pow( 10.0 , get_bad_Q(data,lower, upper)) ;
		  improved_upper_limit_plus_1 = pow( 10.0 , get_bad_Q(data,lower, upper) ) ;
		}
		
		cout << "data.MChi = " << data.MChi << endl;
		cout << "previous upper limt = " << upper << " improved upper limit = " << improved_upper_limit << endl;
		
	  F.function = &get_min_Q;
	  
	  double minQ_iterative = find_min_Q(F, lower, improved_upper_limit, data, true);
		
		cout << "minimum delta M " << get_min_Q(minQ_iterative,&data) << " at Q = " << minQ_iterative << endl;
	
		decays_iterative << " " << decays.calc_lifetime(get_min_Q(minQ_iterative,&data));
	
		// find the maximum mass splitting
	  F.function = &get_max_Q;
		double maxQ_iterative = find_min_Q(F, lower, improved_upper_limit_plus_1, data, false);
				
		cout << "maximum delta M " << -get_max_Q(maxQ_iterative,&data) << " at Q = " << maxQ_iterative << endl;
	
		decays_iterative << " " << decays.calc_lifetime(-get_max_Q(maxQ_iterative,&data));
	
		// non-iterative calculation
		
		cout << "--- non-iterative --- " << endl;
		
	  F.function = &get_min_Q_low;
	  
	  double minQ_non_iterative = find_min_Q(F, lower, upper, data, false);
		
		cout << "minimum delta M " << get_min_Q_low(minQ_non_iterative,&data) << " at Q = " << minQ_non_iterative << endl;
	
		decays_explicit << " " << decays.calc_lifetime(get_min_Q_low(minQ_non_iterative,&data));
	
		// find the maximum mass splitting
	  F.function = &get_max_Q_low;
	  
	  double maxQ_non_iterative = find_min_Q(F, lower, upper, data, false);
		cout << "maximum delta M " << -get_max_Q_low(maxQ_non_iterative,&data) << " at Q = " << maxQ_non_iterative << endl;

		decays_explicit << " " << decays.calc_lifetime(-get_max_Q_low(maxQ_non_iterative,&data));

		// plot with respect to pole mass, factor of 2 involved in FS model definition
		myfile << data.MChi << " " << get_min_Q(minQ_iterative,&data) << " " << -get_max_Q(maxQ_iterative,&data) << " " << get_min_Q_low(minQ_non_iterative,&data) << " " << -get_max_Q_low(maxQ_non_iterative,&data) << endl;
		
		decays_explicit << endl;
		decays_iterative << endl;
		
	}
	
	myfile.close();
	decays_explicit.close();
	decays_iterative.close();
	
}

void Figures::plot_Q(Data data)
{
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings_Q.txt");
  int pts = 15;
  double M = 1000;
  int status = 0;
  double max_Q = 400; // (GeV)
  double min_Q = 50; // (GeV)
  data.MChi = M;
  data.P = M;
  for (int i = 0; i < pts+1 ; i++)
  {
    data.mt = 163.3;  
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    MSSM mssm_1loop(data);
    mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    double delta_m_1 = mssm_1loop.data.SE_1["F12_g1"] - mssm_1loop.data.SE_1["F11_g1"];
    MSSM mssm_2loop(data);
    mssm_2loop.compute_spectra_MB_2loop();
    mssm_2loop.compute_tsil();
    double delta_m_2 = (mssm_2loop.data.SE_1["F12_g1"] - mssm_2loop.data.SE_1["F11_g1"]) + ( mssm_2loop.data.SE_2["F12_g1"] - mssm_2loop.data.SE_2["F11_g1"]);
    myfile << pow(data.Q,0.5) << " " << delta_m_1 <<  " " << delta_m_2 << endl;
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
void Figures::plot_M(Data data)
{
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 100;
  double n = 0;
  int status = 0;
  double max_M = 4000; // (GeV)
  double min_M = 90; // (GeV)
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  data.mt = 163.3;
  data.Q = pow(163.3,2);
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    data.MChi = pow(10,n);
    data.P = data.MChi;
    MSSM mssm_1loop(data);
    mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    double delta_m_1 = mssm_1loop.data.SE_1["F12_g1"] - mssm_1loop.data.SE_1["F11_g1"];
    MSSM mssm_2loop(data);
    mssm_2loop.compute_spectra_MB_2loop();
    mssm_2loop.compute_tsil();
    double delta_m_2 = (mssm_2loop.data.SE_1["F12_g1"] - mssm_2loop.data.SE_1["F11_g1"]) + ( mssm_2loop.data.SE_2["F12_g1"] - mssm_2loop.data.SE_2["F11_g1"]);
    myfile << data.MChi << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  myfile.close();
}

void Figures::plot_M_flexiblesusy(Data data)
{ 
	data.do_tsil_all = false;
  ofstream mass_splittings_iterative;
  mass_splittings_iterative.open ("mass_splittings/data/mass_splittings_iterative.txt");
 
  ofstream pole_mass_n_iterative;
  pole_mass_n_iterative.open ("mass_splittings/data/pole_mass_n_iterative.txt");
 
	ofstream pole_mass_c_iterative;
  pole_mass_c_iterative.open ("mass_splittings/data/pole_mass_c_iterative.txt");
  
  ofstream mass_splittings_explicit;
  mass_splittings_explicit.open ("mass_splittings/data/mass_splittings_explicit.txt");
 
  ofstream pole_mass_n_explicit;
  pole_mass_n_explicit.open ("mass_splittings/data/pole_mass_n_explicit.txt");
 
	ofstream pole_mass_c_explicit;
  pole_mass_c_explicit.open ("mass_splittings/data/pole_mass_c_explicit.txt");  
 
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 500;
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
			MSSM mssm(data);
			mssm.compute_spectra_flexiblesusy();
			
			MSSM mssm_iterative = mssm;
      mssm_iterative.compute_tsil_iterative();
			
			mass_splittings_iterative << " " << mssm_iterative.data.SE_1["F12_g1"] - mssm_iterative.data.SE_1["F11_g1"];
			
			pole_mass_n_iterative << " " << M + mssm_iterative.data.SE_1["F11_g1"];
			pole_mass_c_iterative << " " << M + mssm_iterative.data.SE_1["F12_g1"];
			
			MSSM mssm_explicit = mssm;
      mssm_explicit.compute_tsil();
			
			mass_splittings_explicit << " " << mssm_explicit.data.SE_1["F12_g1"] - mssm_explicit.data.SE_1["F11_g1"];
			
			pole_mass_n_explicit << " " << M + mssm_explicit.data.SE_1["F11_g1"];
			pole_mass_c_explicit << " " << M + mssm_explicit.data.SE_1["F12_g1"];
			
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


void Figures::plot_M_flexiblesusy_2loop(Data data)
{ 
  ofstream mass_splittings_explicit;
  mass_splittings_explicit.open ("mass_splittings/data/mass_splittings_explicit.txt");
 
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 20;
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
    
    mass_splittings_explicit << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];
			MSSM mssm(data);
			mssm.compute_spectra_flexiblesusy();
			
			MSSM mssm_explicit = mssm;
      mssm_explicit.compute_tsil();
			
			mass_splittings_explicit << " " << mssm_explicit.data.SE_1["F12_g1"]+mssm_explicit.data.SE_2["F12_g1"] - mssm_explicit.data.SE_1["F11_g1"]+mssm_explicit.data.SE_2["F11_g1"];
			
		}
    
		mass_splittings_explicit << endl;
  }
  mass_splittings_explicit.close();
  
}
