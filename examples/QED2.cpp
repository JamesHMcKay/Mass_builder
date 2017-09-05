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


  TSIL_COMPLEXCPP Aa, Af, Bfa, Jfff, Kffa ;
  TSIL_COMPLEXCPP Tfff, Vfffa, Baf ; 
  TSIL_COMPLEXCPP  i;
    double Zeta;
  TSIL_REAL  mf,  mf2 ,  ma,  ma2 ,  null, null2 ;

  TSIL_REAL Pi;
  TSIL_REAL EL, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, Z1, Z2m, Z2z, Z3, Z3m ;

  TSIL_REAL s;
  
  
  int init(Data data)
  {
		dcomp ii=-1;ii=sqrt(ii);i=ii;
    Zeta=PI*PI/6;
    Pi=PI;
    mf = data.mf, mf2 = TSIL_POW(data.mf, 2) ,     ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,     null = data.null , null2 = TSIL_POW(data.null, 2) ;

    EL = data.EL;
    Q2 = data.Q;
    p=data.P;
    
    s = TSIL_POW(p,2);
    
    // evaluate TSIL integrals
    Aa = -i*TSIL_A_ (ma2 , Q2);
	
	  Af = -i*TSIL_A_ (mf2 , Q2);
	
	  Bfa = i*TSIL_B_ (mf2, ma2, s, Q2);
	  Baf = Bfa;
	
	  TSIL_SetParametersST (&bar,mf2, mf2, mf2, Q2);
	  TSIL_Evaluate (&bar, s);
	  Jfff= TSIL_GetFunction (&bar,"Suxv");
	
	  Kffa = TSIL_I2_(mf2, mf2, ma2, Q2);
	
	  TSIL_SetParametersST (&bar,mf2, mf2, mf2, Q2);
	  TSIL_Evaluate (&bar, s);
	  Tfff= -TSIL_GetFunction (&bar,"Txuv");
	
	  TSIL_SetParameters (&bar,ma2, mf2, mf2 , 1.0 , mf2, Q2);
	  TSIL_Evaluate (&bar, s);
	  Vfffa= -TSIL_GetFunction (&bar,"Uzxyv");

    return 0;
  }
  
 
 
 
 // calculate diagram 2 with and without the onshell condition for a 
 // zero photon mass as input
 
TSIL_COMPLEXCPP  onshell(Data data)
{
	init(data);
  TSIL_COMPLEXCPP C0 = (Power(EL,4)*(-6*Power(ma,4)*Power(mf,3) - 57*Power(ma,2)*Power(mf,5) + 30*Power(mf,7) + 5*Power(ma,4)*Power(mf,2)*p - 11*Power(ma,2)*Power(mf,4)*p - 42*Power(mf,6)*p + 12*Power(ma,2)*(Power(mf,3) + Power(ma,2)*p - 3*Power(mf,2)*p)*Ae(ma) - 12*Power(ma,2)*(2*Power(mf,3) + Power(ma,2)*p - 2*Power(mf,2)*p)*Ae(mf) + Complex(0,60)*Power(ma,4)*Power(mf,3)*Baf - Complex(0,24)*Power(ma,2)*Power(mf,5)*Baf - Complex(0,96)*Power(mf,7)*Baf + Complex(0,28)*Power(ma,6)*p*Baf - Complex(0,60)*Power(ma,4)*Power(mf,2)*p*Baf - Complex(0,24)*Power(ma,2)*Power(mf,4)*p*Baf + Complex(0,32)*Power(mf,6)*p*Baf + 24*Power(ma,4)*Power(mf,3)*Be(mf,ma) - 72*Power(ma,2)*Power(mf,5)*Be(mf,ma) + 12*Power(ma,6)*p*Be(mf,ma) - 48*Power(ma,4)*Power(mf,2)*p*Be(mf,ma) + 24*Power(ma,2)*Power(mf,4)*p*Be(mf,ma)))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2)));
  TSIL_COMPLEXCPP CAa =   (4*Power(EL,4)*(Complex(0,-1)*(2*Power(ma,2)*Power(mf,2)*(3*mf - 11*p) + 7*Power(ma,4)*p + 6*Power(mf,4)*(-mf + p)) - 6*(Power(mf,3) + Power(ma,2)*p - Power(mf,2)*p)*Af))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CAf =   (4*Power(EL,4)*(Complex(0,1)*(15*Power(ma,2)*Power(mf,2)*(mf - p) + 7*Power(ma,4)*p + 2*Power(mf,4)*(9*mf + p)) - 6*(Power(mf,3) + Power(ma,2)*p - Power(mf,2)*p)*Aa + 6*Power(ma,2)*(2*Power(mf,3) + Power(ma,2)*p - 2*Power(mf,2)*p)*Baf))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CBfa =   (8*Power(EL,4)*(Complex(0,-1)*(Power(ma,2)*Power(mf,4)*(15*mf - 17*p) + 7*Power(ma,4)*Power(mf,2)*p + 4*Power(mf,6)*(-3*mf + p)) + 3*Power(ma,2)*(2*Power(mf,3) + Power(ma,2)*p - 2*Power(mf,2)*p)*Af))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CJfff =   (-4*Power(EL,4)*(2*Power(ma,2)*Power(mf,2)*(3*mf - p) + 4*Power(mf,4)*(3*mf - p) + 3*Power(ma,4)*p))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CKffa =   (4*Power(EL,4)*(Power(ma,2)*Power(mf,2)*(mf - p) + 2*Power(mf,4)*(mf - p) + Power(ma,4)*p))/(Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CVfffa =   (4*Power(EL,4)*(2*Power(ma,2)*Power(mf,2)*(mf - p) + Power(ma,4)*p - 2*Power(mf,4)*(mf + p)))/(-(Power(ma,2)*Power(mf,2)) + 4*Power(mf,4)) ;
  TSIL_COMPLEXCPP CAaAf =   -(Power(EL,4)*(-24*Power(mf,3) - 24*Power(ma,2)*p + 24*Power(mf,2)*p))/(6.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CAfAa =   -(Power(EL,4)*(-24*Power(mf,3) - 24*Power(ma,2)*p + 24*Power(mf,2)*p))/(6.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CAfAf =   (Power(EL,4)*(48*Power(mf,3) + 24*Power(ma,2)*p))/(3.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CAfBfa =   -(Power(EL,4)*(48*Power(ma,2)*Power(mf,3) + 24*Power(ma,4)*p - 48*Power(ma,2)*Power(mf,2)*p))/(6.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  TSIL_COMPLEXCPP CBfaAf =   -(Power(EL,4)*(48*Power(ma,2)*Power(mf,3) + 24*Power(ma,4)*p - 48*Power(ma,2)*Power(mf,2)*p))/(6.*Power(ma,2)*Power(mf,2)*(Power(ma,2) - 4*Power(mf,2))) ;
  return  + C0  + Aa * CAa + Af * CAf + Bfa * CBfa + Jfff * CJfff + Kffa * CKffa + Vfffa * CVfffa + Aa * Af * CAaAf + Af * Aa * CAfAa + Af * Af * CAfAf + Af * Bfa * CAfBfa + Bfa * Af * CBfaAf;
}
  
     

TSIL_COMPLEXCPP  offshell(Data data)
{
	init(data);
  TSIL_COMPLEXCPP C0 = (Power(EL,4)*(216*Power(ma,6)*Power(mf,6)*p - 216*Power(ma,4)*Power(mf,8)*p + 432*Power(ma,2)*Power(mf,10)*p - 432*Power(mf,12)*p + Complex(0,56)*Power(ma,12)*p*Baf - Complex(0,228)*Power(ma,10)*Power(mf,2)*p*Baf + Complex(0,108)*Power(ma,8)*Power(mf,4)*p*Baf + Complex(0,164)*Power(ma,6)*Power(mf,6)*p*Baf + Complex(0,12)*Power(ma,4)*Power(mf,8)*p*Baf - Complex(0,144)*Power(ma,2)*Power(mf,10)*p*Baf + Complex(0,32)*Power(mf,12)*p*Baf + 24*Power(ma,12)*p*Be(mf,ma) - 144*Power(ma,10)*Power(mf,2)*p*Be(mf,ma) + 216*Power(ma,8)*Power(mf,4)*p*Be(mf,ma) - 96*Power(ma,6)*Power(mf,6)*p*Be(mf,ma) - 12*Power(ma,10)*mf*(-Power(p,2)) - 120*Power(ma,8)*Power(mf,3)*(-Power(p,2)) + 372*Power(ma,6)*Power(mf,5)*(-Power(p,2)) + 384*Power(ma,4)*Power(mf,7)*(-Power(p,2)) + 672*Power(ma,2)*Power(mf,9)*(-Power(p,2)) + 10*Power(ma,10)*p*(-Power(p,2)) - 47*Power(ma,8)*Power(mf,2)*p*(-Power(p,2)) - 360*Power(ma,6)*Power(mf,4)*p*(-Power(p,2)) + 503*Power(ma,4)*Power(mf,6)*p*(-Power(p,2)) + 182*Power(ma,2)*Power(mf,8)*p*(-Power(p,2)) + 1368*Power(mf,10)*p*(-Power(p,2)) + Complex(0,120)*Power(ma,10)*mf*Baf*(-Power(p,2)) - Complex(0,360)*Power(ma,8)*Power(mf,3)*Baf*(-Power(p,2)) - Complex(0,144)*Power(ma,6)*Power(mf,5)*Baf*(-Power(p,2)) + Complex(0,384)*Power(ma,4)*Power(mf,7)*Baf*(-Power(p,2)) - Complex(0,116)*Power(ma,10)*p*Baf*(-Power(p,2)) + Complex(0,264)*Power(ma,8)*Power(mf,2)*p*Baf*(-Power(p,2)) + Complex(0,188)*Power(ma,6)*Power(mf,4)*p*Baf*(-Power(p,2)) - Complex(0,256)*Power(ma,4)*Power(mf,6)*p*Baf*(-Power(p,2)) + Complex(0,192)*Power(ma,2)*Power(mf,8)*p*Baf*(-Power(p,2)) - Complex(0,128)*Power(mf,10)*p*Baf*(-Power(p,2)) + 48*Power(ma,10)*mf*Be(mf,ma)*(-Power(p,2)) - 264*Power(ma,8)*Power(mf,3)*Be(mf,ma)*(-Power(p,2)) + 312*Power(ma,6)*Power(mf,5)*Be(mf,ma)*(-Power(p,2)) - 96*Power(ma,4)*Power(mf,7)*Be(mf,ma)*(-Power(p,2)) - 48*Power(ma,10)*p*Be(mf,ma)*(-Power(p,2)) + 192*Power(ma,8)*Power(mf,2)*p*Be(mf,ma)*(-Power(p,2)) + 54*Power(ma,8)*mf*Power((-Power(p,2)),2) + 186*Power(ma,6)*Power(mf,3)*Power((-Power(p,2)),2) - 684*Power(ma,4)*Power(mf,5)*Power((-Power(p,2)),2) - 816*Power(ma,2)*Power(mf,7)*Power((-Power(p,2)),2) - 15*Power(ma,8)*p*Power((-Power(p,2)),2) + 148*Power(ma,6)*Power(mf,2)*p*Power((-Power(p,2)),2) + 155*Power(ma,4)*Power(mf,4)*p*Power((-Power(p,2)),2) - 738*Power(ma,2)*Power(mf,6)*p*Power((-Power(p,2)),2) - 1512*Power(mf,8)*p*Power((-Power(p,2)),2) - Complex(0,168)*Power(ma,8)*mf*Baf*Power((-Power(p,2)),2) + Complex(0,96)*Power(ma,6)*Power(mf,3)*Baf*Power((-Power(p,2)),2) + Complex(0,384)*Power(ma,4)*Power(mf,5)*Baf*Power((-Power(p,2)),2) + Complex(0,60)*Power(ma,8)*p*Baf*Power((-Power(p,2)),2) - Complex(0,100)*Power(ma,6)*Power(mf,2)*p*Baf*Power((-Power(p,2)),2) - Complex(0,40)*Power(ma,4)*Power(mf,4)*p*Baf*Power((-Power(p,2)),2) + Complex(0,32)*Power(ma,2)*Power(mf,6)*p*Baf*Power((-Power(p,2)),2) + Complex(0,192)*Power(mf,8)*p*Baf*Power((-Power(p,2)),2) - 72*Power(ma,8)*mf*Be(mf,ma)*Power((-Power(p,2)),2) + 240*Power(ma,6)*Power(mf,3)*Be(mf,ma)*Power((-Power(p,2)),2) + 192*Power(ma,4)*Power(mf,5)*Be(mf,ma)*Power((-Power(p,2)),2) + 24*Power(ma,8)*p*Be(mf,ma)*Power((-Power(p,2)),2) - 96*Power(ma,6)*Power(mf,2)*p*Be(mf,ma)*Power((-Power(p,2)),2) - 42*Power(ma,6)*mf*Power((-Power(p,2)),3) + 60*Power(ma,4)*Power(mf,3)*Power((-Power(p,2)),3) + 144*Power(ma,2)*Power(mf,5)*Power((-Power(p,2)),3) - 111*Power(ma,4)*Power(mf,2)*p*Power((-Power(p,2)),3) + 114*Power(ma,2)*Power(mf,4)*p*Power((-Power(p,2)),3) + 648*Power(mf,6)*p*Power((-Power(p,2)),3) + Complex(0,48)*Power(ma,6)*mf*Baf*Power((-Power(p,2)),3) + Complex(0,4)*Power(ma,6)*p*Baf*Power((-Power(p,2)),3) + Complex(0,32)*Power(ma,4)*Power(mf,2)*p*Baf*Power((-Power(p,2)),3) - Complex(0,64)*Power(ma,2)*Power(mf,4)*p*Baf*Power((-Power(p,2)),3) - Complex(0,128)*Power(mf,6)*p*Baf*Power((-Power(p,2)),3) + 24*Power(ma,6)*mf*Be(mf,ma)*Power((-Power(p,2)),3) - 96*Power(ma,4)*Power(mf,3)*Be(mf,ma)*Power((-Power(p,2)),3) + 5*Power(ma,4)*p*Power((-Power(p,2)),4) + 10*Power(ma,2)*Power(mf,2)*p*Power((-Power(p,2)),4) - 72*Power(mf,4)*p*Power((-Power(p,2)),4) - Complex(0,4)*Power(ma,4)*p*Baf*Power((-Power(p,2)),4) - Complex(0,16)*Power(ma,2)*Power(mf,2)*p*Baf*Power((-Power(p,2)),4) + Complex(0,32)*Power(mf,4)*p*Baf*Power((-Power(p,2)),4) - 24*Power(ma,6)*(Power(ma,2) - 4*Power(mf,2))*Ae(mf)*((Power(ma,2) - Power(mf,2))*p + (2*mf - p)*(-Power(p,2))) + 24*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Ae(ma)*(Power(ma,2)*(Power(ma,2) - Power(mf,2))*p + (Power(ma,2)*(mf - 2*p) + Power(mf,2)*(mf - p))*(-Power(p,2)) + (-mf + p)*Power((-Power(p,2)),2))))/(6.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*(-Power(p,2))*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*(-Power(p,2)) + Power((-Power(p,2)),2)));
  TSIL_COMPLEXCPP CAa =   (Complex(0,-0.6666666666666666)*Power(EL,4)*(-(p*Power(Power(p,2),3)*(Power(ma,4) + 2*Power(ma,2)*Power(mf,2) - 12*Power(mf,4) - Complex(0,12)*Power(mf,2)*Af)) + (Power(ma,2) - Power(mf,2))*p*(14*Power(ma,8) - 57*Power(ma,6)*Power(mf,2) + 21*Power(ma,4)*Power(mf,4) + 10*Power(ma,2)*Power(mf,6) + 12*Power(mf,8) - Complex(0,12)*(Power(ma,6) - 2*Power(ma,4)*Power(mf,2) - Power(ma,2)*Power(mf,4) - Power(mf,6))*Af) - Power(Power(p,2),2)*(4*Power(ma,6)*(3*mf - 4*p) + 6*Power(ma,2)*Power(mf,4)*(8*mf - p) + 36*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(-60*mf + 49*p) - Complex(0,12)*(-4*Power(ma,2)*Power(mf,3) + Power(ma,4)*(mf - p) - 3*Power(mf,4)*p)*Af) + Power(p,2)*(Power(ma,8)*(12*mf - 29*p) + 6*Power(ma,2)*Power(mf,6)*(8*mf - p) + 36*Power(mf,8)*p + Power(ma,4)*Power(mf,4)*(-12*mf + 13*p) + Power(ma,6)*(-48*Power(mf,3) + 94*Power(mf,2)*p) - Complex(0,12)*(-4*Power(ma,2)*Power(mf,5) + Power(ma,6)*(mf - 2*p) - 3*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(-3*mf + 2*p))*Af)))/(Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CAf =   (2*Power(EL,4)*(-12*Power(mf,2)*p*Power(Power(p,2),4)*Baf + (Power(ma,2) - Power(mf,2))*p*(Complex(0,1)*(14*Power(ma,8) - 57*Power(ma,6)*Power(mf,2) - 33*Power(ma,4)*Power(mf,4) + 10*Power(ma,2)*Power(mf,6) - 96*Power(mf,8)) - 12*(Power(ma,6) - 2*Power(ma,4)*Power(mf,2) - Power(ma,2)*Power(mf,4) - Power(mf,6))*Aa + 12*(Power(ma,8) - 3*Power(ma,6)*Power(mf,2) + Power(ma,4)*Power(mf,4) + Power(mf,8))*Baf) + Power(Power(p,2),3)*(Complex(0,1)*(Power(ma,4) + 26*Power(ma,2)*Power(mf,2) - 96*Power(mf,4))*p + 12*Power(mf,2)*p*Aa + 12*mf*(Power(ma,4) + 4*Power(mf,3)*p + Power(ma,2)*mf*(2*mf + p))*Baf) - Power(Power(p,2),2)*(Complex(0,1)*mf*(6*Power(ma,6) - 288*Power(mf,5)*p - 2*Power(ma,2)*Power(mf,3)*(72*mf + 41*p) + Power(ma,4)*mf*(12*mf + 67*p)) - 12*(-4*Power(ma,2)*Power(mf,3) + Power(ma,4)*(mf - p) - 3*Power(mf,4)*p)*Aa + 12*(Power(ma,6)*(3*mf - p) + 6*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(2*mf + p) + Power(ma,2)*Power(mf,4)*(4*mf + p))*Baf) + Power(p,2)*(Complex(0,1)*(15*Power(ma,8)*(2*mf - p) - 288*Power(mf,8)*p - 2*Power(ma,2)*Power(mf,6)*(72*mf + p) + Power(ma,4)*Power(mf,4)*(-132*mf + 7*p) + 2*Power(ma,6)*Power(mf,2)*(-39*mf + 50*p)) - 12*(-4*Power(ma,2)*Power(mf,5) + Power(ma,6)*(mf - 2*p) - 3*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(-3*mf + 2*p))*Aa + 12*(2*Power(ma,8)*(mf - p) + Power(ma,2)*Power(mf,6)*(2*mf - p) + 4*Power(mf,8)*p + Power(ma,4)*Power(mf,4)*(mf + 2*p) + Power(ma,6)*Power(mf,2)*(-5*mf + 3*p))*Baf)))/(3.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CBfa =   (Complex(0,-1.3333333333333333)*Power(EL,4)*(-(Power(mf,2)*p*Power(Power(p,2),4)*(Power(ma,2) + 2*Power(mf,2) + Complex(0,6)*Af)) + Power(Power(ma,2) - Power(mf,2),2)*p*(14*Power(ma,6)*Power(mf,2) - 33*Power(ma,4)*Power(mf,4) - 15*Power(ma,2)*Power(mf,6) - 2*Power(mf,8) + Complex(0,6)*(Power(ma,6) - 2*Power(ma,4)*Power(mf,2) - Power(ma,2)*Power(mf,4) - Power(mf,6))*Af) + mf*Power(Power(p,2),3)*(mf*(6*Power(ma,2)*Power(mf,2)*(2*mf - p) + 8*Power(mf,4)*p + Power(ma,4)*(12*mf + p)) + Complex(0,6)*(Power(ma,4) + 4*Power(mf,3)*p + Power(ma,2)*mf*(2*mf + p))*Af) - Power(Power(p,2),2)*(Power(mf,2)*(3*Power(ma,6)*(14*mf - 5*p) + 4*Power(ma,2)*Power(mf,4)*(6*mf - p) + 12*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(-24*mf + 7*p)) + Complex(0,6)*(Power(ma,6)*(3*mf - p) + 6*Power(mf,6)*p + Power(ma,4)*Power(mf,2)*(2*mf + p) + Power(ma,2)*Power(mf,4)*(4*mf + p))*Af) + Power(p,2)*(Power(mf,2)*(Power(ma,8)*(30*mf - 29*p) + 3*Power(ma,4)*Power(mf,4)*(20*mf - 7*p) + 8*Power(mf,8)*p + 2*Power(ma,2)*Power(mf,6)*(6*mf + 7*p) + 2*Power(ma,6)*Power(mf,2)*(-51*mf + 32*p)) + Complex(0,6)*(2*Power(ma,8)*(mf - p) + Power(ma,2)*Power(mf,6)*(2*mf - p) + 4*Power(mf,8)*p + Power(ma,4)*Power(mf,4)*(mf + 2*p) + Power(ma,6)*Power(mf,2)*(-5*mf + 3*p))*Af)))/(Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CJfff =   (-4*Power(EL,4)*(3*(Power(ma,8) + Power(ma,6)*Power(mf,2) - Power(ma,4)*Power(mf,4) + Power(ma,2)*Power(mf,6) - 2*Power(mf,8))*p + (Power(ma,6)*(6*mf - 5*p) + 2*Power(ma,4)*Power(mf,2)*(9*mf - p) + 18*Power(mf,6)*p + Power(ma,2)*Power(mf,4)*(12*mf + 7*p))*Power(p,2) - (Power(ma,4)*(6*mf - p) + 18*Power(mf,4)*p + Power(ma,2)*Power(mf,2)*(12*mf + 11*p))*Power(Power(p,2),2) + (Power(ma,2) + 6*Power(mf,2))*p*Power(Power(p,2),3)))/(3.*Power(ma,4)*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CKffa =   (4*Power(EL,4)*((Power(ma,10) - 3*Power(ma,8)*Power(mf,2) - 2*Power(ma,6)*Power(mf,4) + 2*Power(ma,4)*Power(mf,6) + 2*Power(mf,10))*p + (-8*Power(ma,2)*Power(mf,7) - Power(ma,6)*Power(mf,2)*(mf - 3*p) + Power(ma,8)*(mf - 2*p) - 6*Power(mf,8)*p + 2*Power(ma,4)*Power(mf,4)*(-5*mf + 4*p))*Power(p,2) + (8*Power(ma,2)*Power(mf,5) + 2*Power(ma,4)*Power(mf,2)*(mf - p) + 6*Power(mf,6)*p + Power(ma,6)*(-mf + p))*Power(Power(p,2),2) - 2*Power(mf,4)*p*Power(Power(p,2),3)))/(Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CTfff =   (-4*Power(EL,4)*Power(mf,2)*(9*Power(mf,4) - 10*Power(mf,2)*Power(p,2) + Power(Power(p,2),2))*((Power(ma,6) - Power(ma,4)*Power(mf,2) + 2*Power(ma,2)*Power(mf,4) - 2*Power(mf,6))*p + (Power(ma,4)*(2*mf - p) + 4*Power(mf,4)*p + 2*Power(ma,2)*Power(mf,2)*(2*mf + p))*Power(p,2) - 2*Power(mf,2)*p*Power(Power(p,2),2)))/(Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CVfffa =   (-4*Power(EL,4)*(Power(Power(ma,2) - Power(mf,2),2)*(Power(ma,8) - 2*Power(ma,6)*Power(mf,2) - 4*Power(ma,4)*Power(mf,4) - 2*Power(ma,2)*Power(mf,6) - 2*Power(mf,8))*p + (-2*Power(ma,4)*Power(mf,6)*(mf - 2*p) + 2*Power(ma,10)*(mf - p) + 2*Power(ma,2)*Power(mf,8)*(2*mf - p) + 8*Power(mf,10)*p + 3*Power(ma,6)*Power(mf,4)*(mf + 2*p) + Power(ma,8)*Power(mf,2)*(-7*mf + 4*p))*Power(p,2) - (Power(ma,8)*(3*mf - p) + 12*Power(mf,8)*p + 2*Power(ma,6)*Power(mf,2)*(-2*mf + p) + 2*Power(ma,4)*Power(mf,4)*(-2*mf + p) + 2*Power(ma,2)*Power(mf,6)*(4*mf + p))*Power(Power(p,2),2) + mf*(Power(ma,6) - 2*Power(ma,4)*Power(mf,2) + 8*Power(mf,5)*p + 2*Power(ma,2)*Power(mf,3)*(2*mf + p))*Power(Power(p,2),3) - 2*Power(mf,4)*p*Power(Power(p,2),4)))/(Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CAaAf =   -(Power(EL,4)*(-48*Power(ma,8)*p + 144*Power(ma,6)*Power(mf,2)*p - 48*Power(ma,4)*Power(mf,4)*p - 48*Power(mf,8)*p - 48*Power(ma,6)*mf*Power(p,2) + 144*Power(ma,4)*Power(mf,3)*Power(p,2) + 192*Power(ma,2)*Power(mf,5)*Power(p,2) + 96*Power(ma,6)*p*Power(p,2) - 96*Power(ma,4)*Power(mf,2)*p*Power(p,2) + 144*Power(mf,6)*p*Power(p,2) + 48*Power(ma,4)*mf*Power(Power(p,2),2) - 192*Power(ma,2)*Power(mf,3)*Power(Power(p,2),2) - 48*Power(ma,4)*p*Power(Power(p,2),2) - 144*Power(mf,4)*p*Power(Power(p,2),2) + 48*Power(mf,2)*p*Power(Power(p,2),3)))/(12.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CAfAa =   -(Power(EL,4)*(-48*Power(ma,8)*p + 144*Power(ma,6)*Power(mf,2)*p - 48*Power(ma,4)*Power(mf,4)*p - 48*Power(mf,8)*p - 48*Power(ma,6)*mf*Power(p,2) + 144*Power(ma,4)*Power(mf,3)*Power(p,2) + 192*Power(ma,2)*Power(mf,5)*Power(p,2) + 96*Power(ma,6)*p*Power(p,2) - 96*Power(ma,4)*Power(mf,2)*p*Power(p,2) + 144*Power(mf,6)*p*Power(p,2) + 48*Power(ma,4)*mf*Power(Power(p,2),2) - 192*Power(ma,2)*Power(mf,3)*Power(Power(p,2),2) - 48*Power(ma,4)*p*Power(Power(p,2),2) - 144*Power(mf,4)*p*Power(Power(p,2),2) + 48*Power(mf,2)*p*Power(Power(p,2),3)))/(12.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CAfAf =   (Power(EL,4)*(48*Power(ma,8)*p - 144*Power(ma,6)*Power(mf,2)*p + 48*Power(ma,4)*Power(mf,4)*p + 48*Power(mf,8)*p + 96*Power(ma,6)*mf*Power(p,2) - 384*Power(ma,4)*Power(mf,3)*Power(p,2) - 48*Power(ma,6)*p*Power(p,2) - 48*Power(ma,4)*Power(mf,2)*p*Power(p,2) - 192*Power(ma,2)*Power(mf,4)*p*Power(p,2) - 144*Power(mf,6)*p*Power(p,2) + 192*Power(ma,2)*Power(mf,2)*p*Power(Power(p,2),2) + 144*Power(mf,4)*p*Power(Power(p,2),2) - 48*Power(mf,2)*p*Power(Power(p,2),3)))/(6.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CAfBfa =   -(Power(EL,4)*(48*Power(ma,10)*p - 192*Power(ma,8)*Power(mf,2)*p + 192*Power(ma,6)*Power(mf,4)*p - 48*Power(ma,4)*Power(mf,6)*p + 48*Power(ma,2)*Power(mf,8)*p - 48*Power(mf,10)*p + 96*Power(ma,8)*mf*Power(p,2) - 240*Power(ma,6)*Power(mf,3)*Power(p,2) + 48*Power(ma,4)*Power(mf,5)*Power(p,2) + 96*Power(ma,2)*Power(mf,7)*Power(p,2) - 96*Power(ma,8)*p*Power(p,2) + 144*Power(ma,6)*Power(mf,2)*p*Power(p,2) + 96*Power(ma,4)*Power(mf,4)*p*Power(p,2) - 48*Power(ma,2)*Power(mf,6)*p*Power(p,2) + 192*Power(mf,8)*p*Power(p,2) - 144*Power(ma,6)*mf*Power(Power(p,2),2) - 96*Power(ma,4)*Power(mf,3)*Power(Power(p,2),2) - 192*Power(ma,2)*Power(mf,5)*Power(Power(p,2),2) + 48*Power(ma,6)*p*Power(Power(p,2),2) - 48*Power(ma,4)*Power(mf,2)*p*Power(Power(p,2),2) - 48*Power(ma,2)*Power(mf,4)*p*Power(Power(p,2),2) - 288*Power(mf,6)*p*Power(Power(p,2),2) + 48*Power(ma,4)*mf*Power(Power(p,2),3) + 96*Power(ma,2)*Power(mf,3)*Power(Power(p,2),3) + 48*Power(ma,2)*Power(mf,2)*p*Power(Power(p,2),3) + 192*Power(mf,4)*p*Power(Power(p,2),3) - 48*Power(mf,2)*p*Power(Power(p,2),4)))/(12.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  TSIL_COMPLEXCPP CBfaAf =   -(Power(EL,4)*(48*Power(ma,10)*p - 192*Power(ma,8)*Power(mf,2)*p + 192*Power(ma,6)*Power(mf,4)*p - 48*Power(ma,4)*Power(mf,6)*p + 48*Power(ma,2)*Power(mf,8)*p - 48*Power(mf,10)*p + 96*Power(ma,8)*mf*Power(p,2) - 240*Power(ma,6)*Power(mf,3)*Power(p,2) + 48*Power(ma,4)*Power(mf,5)*Power(p,2) + 96*Power(ma,2)*Power(mf,7)*Power(p,2) - 96*Power(ma,8)*p*Power(p,2) + 144*Power(ma,6)*Power(mf,2)*p*Power(p,2) + 96*Power(ma,4)*Power(mf,4)*p*Power(p,2) - 48*Power(ma,2)*Power(mf,6)*p*Power(p,2) + 192*Power(mf,8)*p*Power(p,2) - 144*Power(ma,6)*mf*Power(Power(p,2),2) - 96*Power(ma,4)*Power(mf,3)*Power(Power(p,2),2) - 192*Power(ma,2)*Power(mf,5)*Power(Power(p,2),2) + 48*Power(ma,6)*p*Power(Power(p,2),2) - 48*Power(ma,4)*Power(mf,2)*p*Power(Power(p,2),2) - 48*Power(ma,2)*Power(mf,4)*p*Power(Power(p,2),2) - 288*Power(mf,6)*p*Power(Power(p,2),2) + 48*Power(ma,4)*mf*Power(Power(p,2),3) + 96*Power(ma,2)*Power(mf,3)*Power(Power(p,2),3) + 48*Power(ma,2)*Power(mf,2)*p*Power(Power(p,2),3) + 192*Power(mf,4)*p*Power(Power(p,2),3) - 48*Power(mf,2)*p*Power(Power(p,2),4)))/(12.*Power(ma,4)*(Power(ma,2) - 4*Power(mf,2))*Power(p,2)*(Power(Power(ma,2) - Power(mf,2),2) - 2*(Power(ma,2) + Power(mf,2))*Power(p,2) + Power(Power(p,2),2))) ;
  return  + C0  + Aa * CAa + Af * CAf + Bfa * CBfa + Jfff * CJfff + Kffa * CKffa + Tfff * CTfff + Vfffa * CVfffa + Aa * Af * CAaAf + Af * Aa * CAfAa + Af * Af * CAfAf + Af * Bfa * CAfBfa + Bfa * Af * CBfaAf;
}
     
 
 
 
}



int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Data data(options);
  
  ofstream myfile;
  myfile.open ("examples/QED_test.txt");
  int pts = 10;
  double n = 0;
  int status = 0;
  
  double max_M = 2.0; // (GeV)
  double min_M = 0.1; // (GeV)
  
  for (int i = 0; i < pts+1 ; i++)
  {
    data.P = (i)*(max_M - min_M)/pts + min_M;
    data.mf = 1.0;
	  
	  myfile << data.P << " " << real(extra_TSIL_interface::offshell(data));
	  
	  myfile << " "  << real(extra_TSIL_interface::onshell(data)) << endl;
	}
  
  myfile.close();
  
  return 0;
}
