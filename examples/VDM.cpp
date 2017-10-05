/*
 Mass Builder 
 
 James McKay
 May 2017
 
 -- VDM.cpp --
 
 determine mass splittings for vector dark matter triplet model
 
 requires an input list flag at runtime: ./VDM -i models/VDM/input.txt
 */

#include "data.hpp"
#include "self_energy.hpp"
  
#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

double Pi;

using namespace std;


double iterative_mass_V5(Data data)
{
  
  double M = data.MVp;
  
  double diff = 1, precision = 0.0001;
  int iteration = 0, max_iterations = 500;
  
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    double self_energy = data.SE_1["V5"];
    double mass_sqr = M*M + self_energy;
    
    if (mass_sqr < 0){cout << " mass_sqr < 0 " << endl;}
      
    double mass = pow(abs(mass_sqr),0.5);
      
    diff = abs(mass - data.P);

    data.P = mass;
    //data.Q = mass;
    
    iteration++;
  }
  while (diff > precision  && iteration < max_iterations);
  
  if (iteration == max_iterations)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  return data.P;
}

double iterative_mass_V6(Data data)
{
  
  double M = data.MV0;
  
  double diff = 1, precision = 0.0001;
  int iteration = 0, max_iterations = 500;
  
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    double self_energy = data.SE_1["V6"];
    double mass_sqr = M*M + self_energy;
    
    if (mass_sqr < 0){cout << " mass_sqr < 0 " << endl;}
      
    double mass = pow(abs(mass_sqr),0.5);
      
    diff = abs(mass - data.P);

    data.P = mass;
    //data.Q = mass;
    
    iteration++;
  }
  while (diff > precision  && iteration < max_iterations);
  
  if (iteration == max_iterations)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  return data.P;
}


double pole_mass_V5(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  double Mp = pow(abs(data.MVp*data.MVp + (data.SE_1["V5"])),0.5);
  
  if (data.MVp*data.MVp + (data.SE_1["V5"]) < 0)
  {
    cout << "explicit pole mass less than zero!" << endl;
  }
  
  
  return Mp;
}

double pole_mass_V6(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  double Mp = pow(abs(data.MV0*data.MV0 + (data.SE_1["V6"])),0.5);
  
  if (data.MVp*data.MVp + (data.SE_1["V6"]) < 0)
  {
    cout << "explicit pole mass less than zero!" << endl;
  }
  return Mp;
}

// determine MSbar parameters
void compute_spectra(Data &data)
{
  
  double alpha = pow(data.EE,2) / (4.*Pi);
  
  cout << "alpha (mz) = " << alpha ;
  double A = (7./(30.*Pi)); // SM
  double alpha_mz = alpha;
  double mu0 = data.MZ;
  double mu = data.Q;
  
  
  alpha = pow(   1.0/alpha_mz -  A * log( mu/mu0) , -1);
  
  cout << ", alpha(" << data.Q << ") = " << alpha << endl;
  
  data.EE = pow( (4.*Pi) * alpha , 0.5) ;
  
}




int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy se;
  Data data(options);
  ofstream deltam;
  deltam.open ("models/VDM/output/mass_splittings.txt");
  
    // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  double EE_mz = data.EE;
  Pi = PI;
  // number of points to plot
  int pts = 100;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    data.MVp = M;
    data.MV0 = M;
    data.P = M;
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    deltam << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];
			compute_spectra(data);
	    deltam << " " << pole_mass_V5(data) - pole_mass_V6(data);
	    data.EE = EE_mz;
		}
		
		deltam << endl;
  }
  
  
  cout << "VDM mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_VDM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  
	deltam.close();
  
  return 0;
}
