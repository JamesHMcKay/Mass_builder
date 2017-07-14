/*
 Mass Builder 
 
 James McKay
 July 2017
 
 -- QED.cpp --
 
 compute two-loop self energy for massless QED and compare to analytically obtained results
 
 */

#include "data.hpp"
#include "self_energy.hpp"

using namespace std;


double get_fermion_1loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  
  double self_energy =  pow(EL,2)*p  * ( 1 - log(pow(p,2)/Q2)) / (16. * pow(Pi,2) ) ;
  
  return self_energy;
}

std::complex<double> get_fermion_2loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  double Zeta = pow(Pi,2)/6.;
  
  double C = pow(EL,4) * p / pow(16. * pow(Pi,2) ,2 );
  
  std::complex<double> i;
  dcomp ii=-1;ii=sqrt(ii);i=ii;
  
  std::complex<double> s(-pow(p,2)/Q2 ,-0.0);
  
  std::complex<double> self_energy_1 = C * ( 5.25 - Zeta - 3.*log(s) + 2.*pow(log(s),2) );
  
  std::complex<double> self_energy_2 = C * (-31.+4.*Zeta -4.*log(s)*(-5. + 2.*log(s)))/8.;
  
  std::complex<double> self_energy_3 = C * ( -3.5 + 2.*log(s) );
  
  std::complex<double> ct_1 = C * (4. - Zeta - 2.*log(s) + pow(log(s),2) ) / 2.;
  
  std::complex<double> ct_34 = C * ( -4. + Zeta + 2.*log(s)-pow(log(s),2) );
  
  return self_energy_1 + self_energy_2 + self_energy_3 + ct_34 + ct_1;
}

double get_photon_1loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  std::complex<double> i;
  dcomp ii=-1;ii=sqrt(ii);i=ii;
  
  std::complex<double> s(-pow(p,2)/Q2 ,-0.0);
  
  std::complex<double> self_energy = 0;
  
  for (int i = 0; i<1; i++)
  {
     self_energy += pow(EL,2)*p*p  * ( 20./9.- (4./3.)*log(s)) / (16. * pow(Pi,2)  ) ;
  }
  
  return real(self_energy);
}

int main()
{
  Data data;
  data.P = 10.;
  data.Q = 100.;
  data.EL = 0.1;
  data.mf = 0.00001;
  data.ma = 0.00001;
  cout.precision(17);
  
  Self_energy self_energy;
  
  self_energy.run_tsil(data);
  
  double Fermion1loopAnalytic = get_fermion_1loop(data);
  
  std::complex<double> Fermion2loopAnalytic = get_fermion_2loop(data);
  
  cout << "Self energy computed analytically = " << -real(Fermion1loopAnalytic + Fermion2loopAnalytic) << endl;
  
  cout << "Self energy from Mass Builder = " << data.SE_1["F02_g1"] + data.SE_2["F02_g1"] << endl;
  
  cout << "Difference = " << real(Fermion1loopAnalytic + Fermion2loopAnalytic) + data.SE_1["F02_g1"] + data.SE_2["F02_g1"] << endl;
  
  return 0;
}
