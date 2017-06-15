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

double get_photon_1loop(Data data)
{
  double EL = data.EL;
  double p = data.P;
  double Q2 = data.Q;
  double Pi = 4.0L*atan(1.0L);
  
  double self_energy = 0;
  
  for (int i = 0; i<3; i++)
  {
     self_energy += pow(EL,2)*p  * ( 20./9.- (4./3.)*log(pow(p,2)/Q2)) / (16. * pow(Pi,2)  ) ;
  }
  
  return self_energy;
}

int main()
{

  Data data;
  data.P = 10.;
  data.Q = 100.;
  data.EL = 0.1;
  data.MM = 0.;
  data.ME = 0.;
  data.ML = 0.;
  data.MA = 0.;
  
  Self_energy self_energy;
  
  self_energy.run_tsil(data);
  
  double Fermion1loopMB = data.SE_1["F02_g1"];

  double Fermion1loopAnalytic = get_fermion_1loop(data);
  
  
  double Photon1loopMB = data.SE_1["V1"];

  double Photon1loopAnalytic = get_photon_1loop(data);
  
  cout << "Mass Builder self energy = " << Fermion1loopMB << endl;
  cout << "Analytic result = " << Fermion1loopAnalytic << endl;
  
  cout << "Mass Builder self energy = " << Photon1loopMB << endl;
  cout << "Analytic result = " << Photon1loopAnalytic << endl;

  
  return 0;
}
