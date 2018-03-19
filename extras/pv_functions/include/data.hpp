// structure to hold all data and constants
// convention for masses is M_x (<upper case>_<lower case>) in all parts of the program


#ifndef DATA_H
#define DATA_H

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>


using namespace std;


struct Data
{
  public:
  
  
  // default SM parameters
  // masses
  long double M_h=125.66L;
  long double M_z=91.1876L;
  long double M_w=80.385L;
  long double M_t=173.15L;
  long double M_pi=0.139570L; // pion mass
  long double M_u= 1.05658372E-01; // muon mass
  long double M_e=5.10998928E-04;
  long double M_a=0.0L;
  //double Lam_h;
  //double mu_h;
  
  // other constants/parameters
  long double Pi=3.1415926535897932384626433833;
  long double v0=246L;
  long double alpha_s=0.1184;
  long double M_pl=1.2e19;
  

  
  // BSM parameters
  
  long double Lam_hs=0.1;
  long double Lam_s=0;
  
  long double M_chi=1000; // MDM degenerate mass
  
  long double M_chi_pole; // place holder for dm mass during rd calculations
  
  long double M_chin,M_chip,M_chipp; // place holder for neutral, charged and doubly charged components of Chi
  
  
 // double cw2=0.767595;//0.77; // FlexibleSUSY value
//  double sw2=1.0-cw2;
  
  long double sw2=(1-pow(M_w/M_z,2));//0.23126;
  long double cw2=1.0-sw2;
  // SM couplings
  
  long double alpha=1.0/127.93; // at M_z
  
  long double alpha_2=alpha/sw2;
  
  long double g1=0.65;
  
  
  long double g2=pow(alpha_2*4*Pi,0.5);
  
  
  
  
  long double G_F=1.16637876E-05;

  long double f_pi=0.13041;
  long double mz,mw,Pi2=pow(Pi,2);
  long double V_ud=0.97425;
  
  long double Q=100;  // renormalisation scale
  
  int components=5; // number of multiplet components, used in some generalised decay calculations
  int method=1; // method used to calculate decay rates
  long double tol=0.000000001L; //tolerance used when calculating pole masses
  
  
  
  //
  
  long double particle_1_SE;
  long double particle_2_SE;
  
  // constructor
  Data (){};
  
  Data(int argc, char* argv[]) {
  long double param [10];
  std::string name [10]; int i=0;

  if (argc==1)
  {
  cout << "Please enter a data file in the format ./main input.txt, using default values " << endl;
  }
  else
  {
  std::ifstream input(argv[1]);
  std::string line;
  while(getline(input, line)) {
        if (!line.length() || line[0] == '#')
           continue;
        std::istringstream iss(line);
        iss>> name[i] >> param[i];
    
        i=i+1;
     }
  }
  for (int n=0;n<i+1;n++)
  {
  if (name[n]=="lambda_hs")
  {
  Lam_hs=param[n];
  }
  if (name[n]=="M_t")
  {
  M_t=param[n];
  }
  if (name[n]=="M_h")
  {
  M_h=param[n];
  }
  if (name[n]=="M_w")
  {
  M_w=param[n];
  }
  if (name[n]=="M_z")
  {
  M_z=param[n];
  }
  if (name[n]=="M_chi")
  {
  M_chi=param[n];
  }
  if (name[n]=="M_chi_pole")
  {
  M_chi_pole=param[n];
  }
  if (name[n]=="g1")
  {
  g1=param[n];
  }
  if (name[n]=="g2")
  {
  g2=param[n];
  }
  
  if (name[n]=="method")
  {
  method=param[n];
  }
  
  if (name[n]=="components")
  {
  components=param[n];
  }
  
  if (name[n]=="tol")
  {
  tol=param[n];
  }
  
  if (name[n]=="Q")
  {
  Q=param[n];
  }
  
  if (name[n]=="M_a")
  {
  M_a=param[n];
  }
  
  }
  }
  
  
  long double get_deltaM(){return M_chip-M_chin;}
  
};
#endif



