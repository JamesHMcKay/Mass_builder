#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include "mdm.hpp"
#include "pv.hpp"

using namespace std;

dcomp MDM::Sigma_Kp(long double p)
{
  
  B_0 B0(Q);
  B_1 B1(Q);
  
  dcomp s = data.sw2 * B1(p,M,data.M_a) + data.cw2 * B1(p,M,mz) + B1(p,M,mw)+1.0L;
  
  //dcomp s = data.sw2 * B1(p,M,data.M_a) + data.cw2 * B1(p,M,mz) + 5.0 * B1(p,M,mw);
  
  return s*(-2.0L*pow(g,2))/(16.0L*data.Pi2);
  
}


dcomp MDM::Sigma_Kpp(long double p)
{
  B_0 B0(Q);
  B_1 B1(Q);
  
  dcomp s = 2.0L*data.sw2 * B1(p,M,data.M_a) + 2*data.cw2 * B1(p,M,mz) + B1(p,M,mw);
  
  return s*(-4.0L*pow(g,2))/(16.0L*data.Pi2);
}


dcomp MDM::Sigma_Mp(long double p)
{
  B_0 B0(Q);
  B_1 B1(Q);
  
  dcomp s = data.sw2 * B0(p,M,data.M_a) + data.cw2 * B0(p,M,mz) + B0(p,M,mw)-1.0L;
  
  return s*(-4.0*pow(g,2)*M)/(16.0L*data.Pi2);
}


dcomp MDM::Sigma_Mpp(long double p)
{
  B_0 B0(Q);
  B_1 B1(Q);
  
  //dcomp s = sw2 * B0(p,M,0.0) + cw2 * B0(p,M,mz) + B0(p,M,mw)-1.0;
  dcomp s = 4.0L*data.sw2 * B0(p,M,data.M_a) + 4.0L*data.cw2 * B0(p,M,mz) +  2.0L*B0(p,M,mw)-3.0L*flag;
  return s*(-4.0L*pow(g,2)*M)/(16.0*data.Pi2);
  
}


dcomp MDM::Sigma_K0(long double p)
{
  B_0 B0(Q);
  B_1 B1(Q);
  
  return -pow(g,2.0L) * ( 4.0L * B1(p,M,mw) + 2.0L )   /(16.0L*data.Pi2);
  
}

dcomp MDM::Sigma_M0(long double p)
{
  B_0 B0(Q);
  B_1 B1(Q);
  
  return -(pow(g,2.0L)*M) * ( 8.0L * B0(p,M,mw) - 4.0L )/(16.0L*data.Pi2);
}


long double MDM::calculate_pole_mass_n()
{
  long double prev=-1.0L;
  long double MFn_tree=M,MFn_pole=M;
  
  while (abs(MFn_pole-prev)>data.tol)
  {
    prev=MFn_pole;
    MFn_pole=real((MFn_tree-Sigma_M0(MFn_pole)) / ( 1.0L + Sigma_K0(MFn_pole)));
  }
  return MFn_pole;
}


long double MDM::calculate_pole_mass_c()
{
  long double prev=-1.0L;
  long double MFc_tree=M, MFc_pole=M;
  
  while (abs(MFc_pole-prev)>data.tol)
  {
    prev=MFc_pole;
    MFc_pole=real((MFc_tree-Sigma_Mp(MFc_pole)) / ( 1.0L + Sigma_Kp(MFc_pole)));
  }
  
  return MFc_pole;
}


long double MDM::calculate_pole_mass_cc()
{
  long double prev=-1;
  long double MFcc_tree=M, MFcc_pole=M;
  //for (int i=1; i<2; i++)
  while (abs(MFcc_pole-prev)>data.tol)
  {
    prev=MFcc_pole;
    MFcc_pole=real((MFcc_tree-Sigma_Mpp(MFcc_pole)) / ( 1.0L + Sigma_Kpp(MFcc_pole)));
  
  }

  return MFcc_pole;
}


long double MDM::calculate_mass_splitting()
{
  return calculate_pole_mass_c()-calculate_pole_mass_n();
}

long double MDM::calculate_mass_splitting2_simple()
{
  return calculate_pole_mass_cc_simple()-calculate_pole_mass_n_simple();
}

long double MDM::calculate_mass_splitting2()
{
  return calculate_pole_mass_cc()-calculate_pole_mass_n();
}

long double MDM::calculate_pole_mass_n_simple()
{
  return M-real( Sigma_M0(M) + M * Sigma_K0(M));
}

long double MDM::calculate_pole_mass_c_simple()
{
  long double MFc_pole= M-real( Sigma_Mp(M) + M * Sigma_Kp(M));
  return MFc_pole;
}

long double MDM::calculate_pole_mass_cc_simple()
{
  long double MFcc_pole= M-real( Sigma_Mpp(M) + M * Sigma_Kpp(M));
  return MFcc_pole;
}


// calculate derivative from analytical form of B0, using simple numerical derivative
long double MDM::Sigma_M0_der(long double p){
  long double h=0.1;
  long double result;
  
  result = real(-Sigma_M0(p+2.0*h) + 8.0L * Sigma_M0(p+h) - 8.0L * Sigma_M0(p-h) + Sigma_M0(p-2.0L*h))/(12.0L*h);
  
  //result = real(Sigma_M0(p+h)-Sigma_M0(p-h))/(2.0*h);
  
  result =  result * 0.5L * 1.0L/p; // account for the fact we take deriviative wrt p^2
  
  return result;
}


// calculate derivative from analytical form of B0, using simple numerical derivative
long double MDM::Sigma_K0_der(long double p)
{
  long double h=0.1;
  long double result;
  
  result = real(-Sigma_K0(p+2.0L*h) + 8.0L * Sigma_K0(p+h) - 8.0L * Sigma_K0(p-h) + Sigma_K0(p-2.0L*h))/(12.0L*h);
  
  //result = real(Sigma_M0(p+h)-Sigma_M0(p-h))/(2.0*h);
  
  result =  result * 0.5L * 1.0L/p; // account for the fact we take deriviative wrt p^2
  
  return result;
}




long double MDM::set_pole_mass(int use_simple)
{
  long double diff=1.0L,M_pole=0;
  long double Mp=data.M_chi_pole;
  if (Mp==0)
  {
    cout << " Enter a non-zero value for M_chi_pole in the input file to use this function " << endl;
    return M;
  }
  else{
    while (abs(diff)>0.001)
    {
      if (use_simple==1)
      {
        M_pole=calculate_pole_mass_n_simple();
      }
      else
      {
        M_pole=calculate_pole_mass_n();
      }
      diff = M_pole-Mp;
      M=M-diff;
    }
    data.M_chi=M;
    return M;
  }
}


