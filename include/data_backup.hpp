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
  
  /*  PLANS
  
  There is lots that can be done here, mainly through abstration of particle names and
  how this will interact with the generator codes.
  
  Make all particle arbitrary, and refer to them with strings, or tags
  these should correspond to the same tags as used in the generator
  
  for example, we have 10 free particles, p1 - p10, these are always arbitrary
  and can get assigned on the creation of the code.  Need to get these out of the
  FeynArts model file perhaps.  So every entry in this file needs to have a tag
  and a mass value.
  
  Need to have a clear system for defining which mass is being updated by
  the self energy calculation as well
  
  
  */
  
  
  
  // default SM parameters
  // masses
  double M_h=125.66;
  double M_z=91.1876;
  double M_w=80.385;
  double M_t=173.15;
  double M_pi=0.139570; // pion mass
  double M_u= 1.05658372E-01; // muon mass
  double M_e=5.10998928E-04;
  double M_a=0.0;
  //double Lam_h;
  //double mu_h;
  
  // other constants/parameters
  double Pi=3.14159265359;
  double v0=246;
  double alpha_s=0.1184;
  double M_pl=1.2e19;
  

  
  // BSM parameters
  
  double M;
  double P;
  
  double Lam_hs=0.1;
  double Lam_s=0;
  
  double M_chi=1000; // MDM degenerate mass
  
  double M_chi_pole; // place holder for dm mass during rd calculations
  
  double M_chin,M_chip,M_chipp; // place holder for neutral, charged and doubly charged components of Chi
  
  
  double tW = acos(M_w/M_z);
  double sw2 = pow(sin(tW),2),S2TW=pow(sin(tW),2), cw2 = pow(cos(tW),2);
  double sw = sin(tW), cw = cos(tW);
  
  double ms;
  
//  double sw2=(1-pow(M_w/M_z,2));//0.23126;
//  double cw2=1.0-sw2;
  // SM couplings
  
  double alpha=1.0/127.93; // at M_z
  
  double alpha_2=alpha/sw2;
  
  double g1=0.65;
  double lambda;
  double g;
  
  
  double g2=pow(alpha_2*4*Pi,0.5);
  
  
  
  
  double G_F=1.16637876E-05;

  double f_pi=0.13041;
  double mz,mw,Pi2=pow(Pi,2);
  double V_ud=0.97425;
  
  double Q=100;  // renormalisation scale
  
  int components=5; // number of multiplet components, used in some generalised decay calculations
  int method=1; // method used to calculate decay rates
  double tol=0.00001; //tolerance used when calculating pole masses
  
  
  
  //
  
  double SE_1;
  double SE_2;
  
  
  
  
  // constructor
  Data (){};
  
  Data(int argc, char* argv[]) {
  double param [99];
  std::string name [99]; int i=0;

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
  
  if (name[n]=="M")
  {
  M=param[n];
  }
  
  if (name[n]=="P")
  {
  P=param[n];
  }
  
  if (name[n]=="lambda")
  {
  lambda=param[n];
  }
  
  if (name[n]=="g")
  {
  g=param[n];
  }
  
  }
  }
  
  
  double get_deltaM(){return M_chip-M_chin;}
  
};
#endif



