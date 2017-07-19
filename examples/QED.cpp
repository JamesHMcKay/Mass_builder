/*
 Mass Builder 
 
 James McKay
 July 2017
 
 -- QED.cpp --
 
 compute two-loop self energy for massless QED and compare to analytically obtained results
 
 */

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

using namespace std;


// extra TSIL interface for manually entered derivatives
// of one-loop self energies
  
#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

namespace extra_TSIL_interface
{
  #include TSIL_PATH
  using namespace std;
  TSIL_DATA bar;

  long double strtold(const char *, char **);
  // define subroutines here
  TSIL_REAL Q2,Q;
  TSIL_REAL p;
  TSIL_COMPLEXCPP Log(TSIL_REAL a){complex<double> s(a/Q2,-0.000);return log(s);}
  TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}
  TSIL_COMPLEXCPP Power(TSIL_COMPLEXCPP a, int b){return pow(a,b);}
  TSIL_REAL Sin(TSIL_REAL a){return sin(a);}
  TSIL_REAL Cos(TSIL_REAL a){return cos(a);}
  TSIL_REAL Dot(TSIL_REAL a, TSIL_REAL b){return a*b;}
  TSIL_REAL Sqrt(TSIL_REAL a){return TSIL_POW(a,0.5);}
  TSIL_REAL Zeta2 = 1.6449340668482;
  TSIL_COMPLEXCPP Ae(TSIL_REAL a) { return TSIL_Aeps_(TSIL_POW(a,2),Q2);}
  TSIL_COMPLEXCPP Be(TSIL_REAL a, TSIL_REAL b) { return TSIL_Beps_(TSIL_POW(a,2),TSIL_POW(b,2), TSIL_POW(p,2), Q2);}
  int          init(Data data);
  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}
  TSIL_COMPLEXCPP operator+(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}
  TSIL_COMPLEXCPP operator+(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}
  TSIL_COMPLEXCPP operator-(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}
  TSIL_COMPLEXCPP operator-(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}
  TSIL_COMPLEXCPP operator/(TSIL_COMPLEXCPP a,double b){TSIL_COMPLEXCPP c=b;return a/c;}
  TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}
 
  TSIL_COMPLEXCPP Aa, Af, Baf, Bff, Faffaf ;
  TSIL_COMPLEXCPP Jaaf, Jfff, Kaff, Taaf, Tfaa ;
  TSIL_COMPLEXCPP Tfff, Vaaff, Vfaaf, Vfffa, Bfa ;

  TSIL_COMPLEXCPP  i;
    double Zeta;
  TSIL_REAL  ma,  ma2 ,  mf,  mf2 ,  null, null2 ;

  TSIL_REAL Pi;
  TSIL_REAL EL, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, Z1, Z2m, Z2z, Z3, Z3m ;
  
  TSIL_COMPLEXCPP dBfa;
  
  
  
  
  int init(Data data)
  {
    i=Power(-1,0.5);
    std::complex<double> s(pow(data.P,2) ,-0.0);
    Pi = 4.0L*atan(1.0L);
    p = data.P;
    Q2 = data.Q;
    
    ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,   mf = data.mf, mf2 = TSIL_POW(data.mf, 2) ,   null = data.null , null2 = TSIL_POW(data.null, 2) ;
    EL = data.EL;
    
    dBfa = i*TSIL_dBds_(mf2,ma2,TSIL_POW(data.P,2),data.Q);
    
    Baf = i*TSIL_B_ (ma2, mf2, TSIL_POW(data.P,2),data.Q);
    Bfa = Baf;
    return 0;
  }
  
  double fermion_1loop_derivative(Data &data)
  {
    init(data);
    
    TSIL_REAL s = pow(data.P,2);
    
    // try multiplying before taking the order epsilon^0 term
    
    std::complex<long double> result = (Power(EL,4)*(Sqrt(s)*(-2*mf*s + Power(s,1.5) - 4*mf*Ae(ma) + 4*mf*Ae(mf) +
          2*(-2*Power(ma,2)*mf + 2*Power(mf,3) - 4*mf*s + Power(s,1.5))*Be(ma,mf)) + Complex(0,2)*s*Aa -
       Complex(0,2)*s*Af + Complex(0,2)*s*(-Power(ma,2) + Power(mf,2) - 3*mf*Sqrt(s) + s)*
        Bfa + (Power(ma,2) - Power(mf,2) + s)*Bfa*
        (Complex(0,-6)*mf*Sqrt(s) + Complex(0,2)*s - Aa + Af +
          (Power(ma,2) - Power(mf,2) + 4*mf*Sqrt(s) - s)*Bfa) -
       2*s*(Complex(0,-2)*(3*(ma - mf)*mf*(ma + mf)*Sqrt(s) - (Power(ma,2) - 9*Power(mf,2))*s - 6*mf*Power(s,1.5) + Power(s,2)) + 
          (Power(ma,2) - Power(mf,2) + 4*mf*Sqrt(s) - s)*(-Aa + Af +
             (Power(ma,2) - Power(mf,2) + 4*mf*Sqrt(s) - s)*Bfa))*dBfa))/
   (512.*Power(Pi,4)*Power(s,2));
    
    std::complex<double> s2(pow(data.P,2) ,-0.0);
    
    
    // for the massless model the equivalent is
    
    //result = (Power(EL,4)*(-1 + log(-s2))*(3 + log(-s2)))/(512.*Power(Pi,4));
    
    return real(result);
  }
}



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
  
  // print out the epsilon term for B(0,0) for comparison
  
  std::complex<double> Be00 = (s/16.) * ( -4.*Zeta+115.+8.*pow(log(s),2)-52.*log(s));
  
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
  data.mf = 10;
  data.ma = 0;
  cout.precision(17);
  
  double Pi = 4.0L*atan(1.0L);
  
  Self_energy self_energy;
  
  self_energy.run_tsil(data);
  
  double Fermion1loopAnalytic = get_fermion_1loop(data);
  
  std::complex<double> Fermion2loopAnalytic = get_fermion_2loop(data);
  
  cout << "Self energy computed analytically = " << -real(Fermion1loopAnalytic + Fermion2loopAnalytic) << endl;
  
  
  // get derivative of one-loop self energy, as computed by MB
  
  double SE_der = extra_TSIL_interface::fermion_1loop_derivative(data);
  
  
  cout << "SE_der = " << SE_der	 << endl;
  
  cout << "Self energy from Mass Builder = " << data.SE_1["F02_g1"] + data.SE_2["F02_g1"] << endl;
  
  cout << "Difference = " << real(Fermion1loopAnalytic + Fermion2loopAnalytic) + data.SE_1["F02_g1"] + data.SE_2["F02_g1"] << endl;
  
  return 0;
}
