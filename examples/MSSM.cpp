/*
 Mass Builder
 
 -- MSSM.cpp --
 
 example application, determing mass splitting between neutral and charged winos in MSSM
 
 requires an input list flag at runtime: ./MSSM -i models/MSSM/input.txt
 */


#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

using namespace std;

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif



namespace extra_TSIL_interface
{
  #include TSIL_PATH
  using namespace std;

	double B0(double p,double M,double m,double Q)
	{
  
   TSIL_REAL Q2 = pow(Q,2);
   TSIL_REAL s = pow(p,2);

   return real(TSIL_B_ (pow(M,2), pow(m,2), s, Q2));

  }
  
}


double iterative_mass(Data data, string particle)
{
  long double M_tree = data.MChi;
  long double M_pole = data.MChi;
  data.P = M_tree;
  long double diff = 1;
  long double precision = 1e-8;
  int iteration = 0;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    if (data.do_tsil_all)
    {
			M_pole = M_tree + data.SE_1[particle] + data.SE_2[particle];
    }
    else
    {
      M_pole = M_tree + data.SE_1[particle];
    }
    
    diff = abs(M_pole - data.P);
    data.P = M_pole;
    
    #ifdef DEBUG
    cout<< "\r" << "M_pole - p = " << diff << " GeV";
    std::cout << std::flush;
		#endif
		if (iteration > 40 && diff > 1e-3)
		{
			iteration = 49 ; 
		}
    
    iteration++;
  } while (diff > precision  && iteration < 50);
  #ifdef DEBUG
  cout<< "\r" << "M_pole - p = " << diff << " GeV";
  cout << "\n";
  #endif
  
  if (iteration == 50)
  {
		#ifdef DEBUG
    cout << "pole mass did not converge" << endl;
    #endif
    return 0;
  }
    
  return M_pole;
}







void plot_M(Data data)
{
  
  ofstream myfile;
  myfile.open ("mass_splittings.txt");
  int pts = 100;
  double n = 0;
  double M = 0;
  int status = 0;
  
  double max_M = 10000; // (GeV)
  double min_M = 1; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  
  data.do_tsil_all = false;
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    M = pow(10,n);
    data.MChi=M;
    data.P = M;
    
    //data.Q = M;
    
    // compute explicit mass splitting

    // compute iterative mass splitting
    
    double delta_m1_iterative_1loop = iterative_mass(data,"F12_g1") - iterative_mass(data,"F11_g1");
    
    double r = iterative_mass(data,"F12_g1")/ M;
    //cout << "M_pole / M = " << iterative_mass(data,"F12_g1")/ M << endl;
    
    data.P = iterative_mass(data,"F12_g1");
    
    Self_energy se;
    se.run_tsil(data);    
    
    
    
    
    double delta_m1_explicit_1loop = (data.SE_1["F12_g1"]) - (data.SE_1["F11_g1"]);
    
    
    myfile << M << " " << delta_m1_explicit_1loop  <<  " " << delta_m1_iterative_1loop << endl;
    //status=(float(i)/pts)*100;
    //cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    //std::cout << std::flush;
  }
  
  //status=100;
  //cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  //cout << "\n";
  

  
  myfile.close();
}




void plot_limits(Data data)
{
  ofstream myfile;
  myfile.open ("examples/limits.txt");
  using namespace extra_TSIL_interface;
  long double Q=163,q,q1,q2,q3,q4,q5;
  double mw = data.mw;
  double mz = data.mz;
  
	double cw =  mw/mz;
	double cw2 =  TSIL_POW(cw,2);
	
	double sw =  TSIL_POW(1.-cw2,0.5);
	double sw2 =  TSIL_POW(sw,2);
	double Pi = PI;
	double e = pow(4*Pi*data.alpha,0.5);
	
	double c = 2.0*pow(e,2)/(16.0*(3.14159L)*sw2);

  long double n,M;
  for (int i=1;i<100;i++)
  {
    n=(float(i)/float(100))*8;  //(float(i)/100)*2000+1;
    M=pow(10,n);
    
    q1=c*real(M*(sw2*B0(M,M,0,Q)+cw2*B0(M,M,mz,Q)-B0(M,M,mw,Q))/(3.14159L));
    q=0.9999;
    q2=0.5*(4.0-1/q-q)*c*real(M*(sw2*B0(M*q,M,0,Q)+cw2*B0(M*q,M,mz,Q)-B0(M*q,M,mw,Q))/(3.14159L));
    q=0.99;
    q3=0.5*(4.0-1/q-q)*c*real(M*(sw2*B0(M*q,M,0,Q)+cw2*B0(M*q,M,mz,Q)-B0(M*q,M,mw,Q))/(3.14159L));
    q=1.001;
    q4=0.5*(4.0-1/q-q)*c*real(M*(sw2*B0(M*q,M,0,Q)+cw2*B0(M*q,M,mz,Q)-B0(M*q,M,mw,Q))/(3.14159L));
    q=1.0001;
    q5=0.5*(4.0-1/q-q)*c*real(M*(sw2*B0(M*q,M,0,Q)+cw2*B0(M*q,M,mz,Q)-B0(M*q,M,mw,Q))/(3.14159L));
    
    myfile << M << " " << q1 << " " << q2 << " " << q3 << " " << q4 << " " << q5 << endl;
    
  }
  myfile.close();
  //system("python ../Figures/limits.py");
}






int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  Data data(options);
  
  plot_M(data);
  plot_limits(data);
  
  cout << "example mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_MSSM.py "<< endl;
  cout << "          python examples/plot_limits.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  
  return 0;
}
