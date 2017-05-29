/*
 Mass Builder 
 
 James McKay
 May 2017
 
 -- VDM.cpp --
 
 determine mass splittings for vector dark matter triplet model
 
 requires an input list flag at runtime: ./VDM -i models/VDM/input.txt
 */

#include "data.hpp"
#include "calc_amplitudes.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"


using namespace std;
using namespace utils;


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
    data.Q = mass;
    
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
    data.Q = mass;
    
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
  
  return Mp;
}



int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy se;
  Data data(options);
  
  ofstream myfile;
  ofstream myfile2;
  myfile.open ("models/VDM/output/mass_splittings.txt");
  myfile2.open ("models/VDM/output/masses.txt");
  int pts = 300;
  double n = 0;
  double M=0;
  int status = 0;
  for (int i = 0; i < pts ; i++)
  {
    n=(float(i)/float(pts))*5;
    M= pow(10,n)+10;
    data.MVp=M;
    data.MV0=M;
    data.P = M;
    data.Q = pow(M,2);
    
    // calculate iterative deltaM with Q^2 = 100
    double delta_m_it= iterative_mass_V5(data) - iterative_mass_V6(data);
    data.MVp=M;
    data.MV0=M;
    data.P = M;
    
    // calculate deltaM with Q^2 = 100
    
    data.Q = 100;
    double delta_m= pole_mass_V5(data) - pole_mass_V6(data);
    
    // calculate deltaM with Q^2 = M^2
    data.Q = pow(M,2);
    double delta_m2= pole_mass_V5(data) - pole_mass_V6(data);
    
    myfile << M << " " << delta_m_it << " " << delta_m << " " << delta_m2 << endl;
    
    
    myfile2 << M << " " << pole_mass_V5(data) << " " << pole_mass_V6(data) << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  
  cout << "VDM mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python src/plot_VDM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  
  myfile.close();
  myfile2.close();
  
  return 0;
}
