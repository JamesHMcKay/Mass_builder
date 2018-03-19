#ifndef MDM_H
#define MDM_H

#include "pv.hpp"
#include "data.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>

using namespace std;
typedef complex<long double> dcomp;

class MDM
{
private:
  long double Q;
  int flag=0;
  long double mz,mw,M,g;
public:
  Data data;
  MDM (){}  // defualt constructor
  MDM(Data _data)
  {
    data=_data;
    mw=data.M_w;
    mz=data.M_z;
    Q=data.Q;
    g=data.g2;
    M=data.M_chi;
    
  } //constructor
  
  dcomp Sigma_Kp(long double p );
  dcomp Sigma_Kpp(long double p );
  dcomp Sigma_Mp(long double p );
  dcomp Sigma_Mpp(long double p );
  
  dcomp Sigma_K0(long double p );
  dcomp Sigma_M0(long double p );
  
  
  dcomp Sigma_c(long double p)
  {
    return p*Sigma_Kp(p)+Sigma_Mp(p);
  }
  
  
  dcomp Sigma_n(long double p)
  {
    return p*Sigma_K0(p)+Sigma_M0(p);
  }
  
  long double Sigma_M0_der(long double p);
  long double Sigma_K0_der(long double p);
  long double Sigma_M0_der_integral(long double p);
  long double calculate_pole_mass_n();
  long double calculate_pole_mass_c();
  long double calculate_pole_mass_cc();
  long double calculate_pole_mass_n_simple();
  long double calculate_pole_mass_c_simple();
  long double calculate_pole_mass_cc_simple();
  
  
  long double calculate_mass_splitting();
  long double calculate_mass_splitting2();
  long double calculate_mass_splitting2_simple();
  long double set_pole_mass(int use_simple=0);
  
};

#endif