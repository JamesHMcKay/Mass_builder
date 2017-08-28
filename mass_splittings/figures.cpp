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
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  //mssm.compute_tsil_iterative();
	
  double deltam =  (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
	cout << "checking deltam = " << deltam << endl;
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
  //mssm.compute_tsil_iterative();
	
  return (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double get_max_Q(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
  MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  //mssm.compute_tsil_iterative();
	
  return -(mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}


double get_max_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
	
	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  //mssm.compute_tsil();
	
  return -(mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double get_min_Q_low(double scale, void *params)
{
	Data data = *(Data*)params;
	
	data.Q = scale;
  
	MSSM mssm(data);
	mssm.compute_spectra_flexiblesusy();
	//mssm.compute_spectra_MB_1loop();
  //mssm.compute_tsil();
	
  return (mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"]);
}

double find_min_Q(gsl_function F, double lower, double upper,double M)
{
				
	int status;
	int iteration = 0, max_iteration = 1000;
	
	Data data;
	data.MChi = M; // (GeV)
	
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
      

void Figures::plot_uncertainties(Data data)
{
	
	ofstream myfile;
	myfile.open ("mass_splittings/uncertainties.txt");
  
  double n;
  // set range of plot
  long double logMax = log10(10000);
  long double logMin = log10(10.0);
  
  double mt = 173.15;

	double improved_upper_limit;
	double improved_upper_limit_plus_1;  
  
  // number of points to plot
  int pts = 100;
  for (int i=1;i<pts+1;i++)
	{
		
    double n = i*(logMax - logMin)/pts + logMin;		
		
		data.MChi = pow(10.0L,n);
		
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
		
		
		if (data.MChi > 1000 )
		{
			improved_upper_limit = pow( 10.0 , get_bad_Q(data,lower, upper)-1.0 ) ;
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
	  
	  double minQ_iterative = find_min_Q(F, lower, improved_upper_limit, data.MChi);
		
		cout << "minimum delta M " << get_min_Q(minQ_iterative,&data) << " at Q = " << minQ_iterative << endl;
	
		// find the maximum mass splitting
	  F.function = &get_max_Q;
		double maxQ_iterative = find_min_Q(F, lower, improved_upper_limit_plus_1, data.MChi);
				
		cout << "maximum delta M " << -get_max_Q(maxQ_iterative,&data) << " at Q = " << maxQ_iterative << endl;
	
		// non-iterative calculation
		
		cout << "--- non-iterative --- " << endl;
		
	  F.function = &get_min_Q_low;
	  
	  double minQ_non_iterative = find_min_Q(F, lower, upper, data.MChi);
		
		cout << "minimum delta M " << get_min_Q_low(minQ_non_iterative,&data) << " at Q = " << minQ_non_iterative << endl;
	
		// find the maximum mass splitting
	  F.function = &get_max_Q_low;
	  
	  double maxQ_non_iterative = find_min_Q(F, lower, upper, data.MChi);
		cout << "maximum delta M " << -get_max_Q_low(maxQ_non_iterative,&data) << " at Q = " << maxQ_non_iterative << endl;

		// plot with respect to pole mass, factor of 2 involved in FS model definition
		myfile << data.MChi << " " << get_min_Q(minQ_iterative,&data) << " " << -get_max_Q(maxQ_iterative,&data) << " " << get_min_Q_low(minQ_non_iterative,&data) << " " << -get_max_Q_low(maxQ_non_iterative,&data) << endl;
		
	}
	
	myfile.close();
	
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
  // plot MDM mass splittings
  
  ofstream myfile_a;
  myfile_a.open ("mass_splittings/mass_splittings_iterative.txt");
  
  //ofstream myfile_b;
  //myfile_b.open ("mass_splittings/mass_splittings_explicit.txt");
  
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 100;
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
    
    myfile_a << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];//pow(Q[i],2);
			MSSM mssm(data);
			mssm.compute_spectra_flexiblesusy();
			//mssm.compute_spectra_MB_1loop();
      //mssm.compute_tsil_iterative();
			
			myfile_a << " " << mssm.data.SE_1["F12_g1"] - mssm.data.SE_1["F11_g1"];
		}
    
    myfile_a << endl;
    
  }
  myfile_a.close();
  //myfile_b.close();
  //system("python ../Figures/mass_splittings.py");
}
