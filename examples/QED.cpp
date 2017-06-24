/*
 Mass Builder 
 
 James McKay
 Sep 2016
 
 -- EW_triplet.cpp --
 
 compute full two-loop self energy including derivative of 1-loop functions
 
 requires an input list flag at runtime: ./EW_triplet -i models/EW_triplet/input.txt
 */

#include "data.hpp"
#include "compute_amp.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"
#include "EW_triplet.hpp"


using namespace std;

using namespace utils;

double get_fermion_1loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  
  double self_energy =  pow(EL,2)*p  * ( 1 - log(pow(p,2)/Q2)) / (16. * pow(Pi,2) ) ;
  
  return self_energy;
}

double get_fermion_2loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  
  double Zeta = pow(Pi,2)/6.;
  
  double x = pow(p,2)/Q2;
  
  double C = pow(EL,4) * p / pow(16. * pow(Pi,2) ,2 );
  
  double self_energy_1 = C * ( 5.25 - Zeta - 3*log(x) + 2*pow(log(x),2) );
  
  double self_energy_2 = C * ( (-31 + 4*Zeta + 20*log(x) - 8*pow(log(x),2))/8. );
  
  double self_energy_3 = C * ( -3.5 + 2*log(x) );
  
  cout << "2-loop fermion diagram 1 = " << self_energy_1 << endl;  // this should match diagram 1
  cout << "2-loop fermion diagram 2 = " << self_energy_3 << endl;  // this matches diagram 2
  cout << "2-loop fermion diagram 5 = " << self_energy_2 << endl;  // this should match diagram 5
  
  // now use results from ArXiv:0508242
  
  //double g1 =
  
  
  
  return self_energy_1 + self_energy_2 + self_energy_3;
}

double get_photon_1loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  
  double self_energy = 0;
  
  for (int i = 0; i<1; i++)
  {
     self_energy += pow(EL,2)*p*p  * ( 20./9.- (4./3.)*log(pow(p,2)/Q2)) / (16. * pow(Pi,2)  ) ;
  }
  
  return self_energy;
}

int main()
{

  Data data;
  data.P = 10.;
  data.Q = 100.;
  data.EL = 0.1;
  data.mf = 0.001;
  data.ma = 0.001;
  
  Self_energy self_energy;
  
  cout << " -----  Mass Builder results ----- " << endl;
  
  self_energy.run_tsil(data);
  
  //double Fermion1loopMB = data.SE_1["F02_g1"];

  double Fermion1loopAnalytic = get_fermion_1loop(data);
  
  //double Photon1loopMB = data.SE_1["V1"];

  double Photon1loopAnalytic = get_photon_1loop(data);
  
  cout << " -----  Analytic results ----- " << endl;
  
  cout << "1-loop fermion diagram 1 = " << Fermion1loopAnalytic << endl;
  
  
  
  // two-loop self energies
  
  //double Fermion2loopMB = data.SE_2["F02_g1"];
  
  get_fermion_2loop(data);
    
  //cout << "Mass Builder two-loop self energy = " << Fermion2loopMB << endl;
  //cout << "Analytic two-loop result = " << Fermion2loopAnalytic << endl;
  
  cout << "1-loop photon diagram 1 = " << Photon1loopAnalytic << endl;
  

  
  return 0;
}
