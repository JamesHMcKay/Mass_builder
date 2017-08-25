/*
 Mass Builder
 
 -- MSSM.cpp --
 
 example application, determing mass splitting between neutral and charged winos in MSSM
 
 requires an input list flag at runtime: ./MSSM -i models/MSSM/input.txt
 */


#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

using namespace std;

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

double Pi = PI;

namespace extra_TSIL_interface
{
#include TSIL_PATH
  using namespace std;
  TSIL_DATA bar;
  
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
  
  
  
  
  std::complex<long double>  i;
  long double  MChi,  MChi2 ,  ma,  ma2 ,mt , mt2,  mw,  mw2 ,  mz,  mz2 ,  mh,  mh2, mf, mf2 ,  null, null2 ;
  
  double Pi;
  long double alpha,CTW, Ca, Cw1, Cw2, Cz1, Cz2, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, MassBuilderJEpsilon, STW, cw, cw2, d1Z, d1m, d2Z, d2m, dMWsq1, dMZsq1, dZAA1, dZAZ1, dZW1, dZZA1, dZZZ1, dg2, e, g, g1, g2, sw, sw2, v ;
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca,Bac,Bcz,Bcw;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,dBac,M;
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az;
  
  
  void DoTSIL(Data data)
  {
    alpha=data.alpha, cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   STW = data.STW,   CTW = data.CTW,   e = data.e,   d1Z = data.d1Z,   d1m = data.d1m,   d2Z = data.d2Z,   d2m = data.d2m,   dZW1 = data.dZW1,   dMWsq1 = data.dMWsq1,   dMZsq1 = data.dMZsq1,   dZAA1 = data.dZAA1,   dZZA1 = data.dZZA1,   dZAZ1 = data.dZAZ1,   dZZZ1 = data.dZZZ1,   dg2 = data.dg2,   Ca = data.Ca,   Cz1 = data.Cz1,   Cz2 = data.Cz2,   Cw1 = data.Cw1,   Cw2 = data.Cw2,   v = data.v,   MassBuilderJEpsilon = data.MassBuilderJEpsilon ;
    
    TSIL_REAL Q2 = data.Q;
    TSIL_REAL s = pow(data.P,2);
    MassBuilderCTM1 = 0;
    MassBuilderCTZ1 = 0;
    MassBuilderCTM2 = 0;
    MassBuilderCTZ1 = 0;
    MChi = data.MChi, MChi2 = TSIL_POW(data.MChi, 2) ,   ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,   mw = data.mw, mw2 = TSIL_POW(data.mw, 2) ,   mz = data.mz, mz2 = TSIL_POW(data.mz, 2) ,   mh = data.mh, mh2 = TSIL_POW(data.mh, 2) ,   null = data.null, null2 = TSIL_POW(data.null, 2) ,   mf = data.mf, mf2 = TSIL_POW(data.mf, 2) ,   mt = data.mt , mt2 = TSIL_POW(data.mt, 2) ;
    
    dcomp ii=-1;ii=sqrt(ii);i=ii;
    Pi=PI;
    cw =  mw/mz;
    cw2 =  TSIL_POW(cw,2);
    sw =  TSIL_POW(1.-cw2,0.5);
    sw2 =  TSIL_POW(sw,2);
    STW =  sw;
    CTW =  cw;
    e =  TSIL_POW(4*Pi*alpha,0.5);
    d1Z =  (2*TSIL_POW(e,2)*(TSIL_POW(-1 + TSIL_POW(sw,2),2) + TSIL_POW(cw,2)*(1 + TSIL_POW(sw,2))))/(TSIL_POW(cw,2)*TSIL_POW(sw,2));
    d1m =  (-8*TSIL_POW(e,2)*MChi*(TSIL_POW(-1 + TSIL_POW(sw,2),2) + TSIL_POW(cw,2)*(1 + TSIL_POW(sw,2))))/(TSIL_POW(cw,2)*TSIL_POW(sw,2));
    d2Z =  (4*TSIL_POW(e,2))/TSIL_POW(sw,2);
    d2m =  (-16*TSIL_POW(e,2)*MChi)/TSIL_POW(sw,2);
    dZW1 =  (TSIL_POW(e,2)*(-9 + 20*TSIL_POW(cw,2) + 20*TSIL_POW(sw,2)))/(6.*TSIL_POW(sw,2));
    dMWsq1 =  -(TSIL_POW(e,2)*(TSIL_POW(cw,4)*(44*TSIL_POW(mw,2) + 6*TSIL_POW(mz,2)) - 6*TSIL_POW(mw,2)*TSIL_POW(sw,4) + TSIL_POW(cw,2)*(6*TSIL_POW(ma,2)*TSIL_POW(sw,2) + TSIL_POW(mw,2)*(-33 + 38*TSIL_POW(sw,2)))))/(6.*TSIL_POW(cw,2)*TSIL_POW(sw,2));
    dMZsq1 =  (TSIL_POW(e,2)*(TSIL_POW(cw,2)*TSIL_POW(mz,2)*(17 - 39*TSIL_POW(cw,4) - 2*(16 + TSIL_POW(cw,2))*TSIL_POW(sw,2) + 17*TSIL_POW(sw,4)) + 12*TSIL_POW(mw,2)*(1 - 2*TSIL_POW(cw,6) + 2*TSIL_POW(cw,2)*TSIL_POW(sw,4))))/(12.*TSIL_POW(cw,4)*TSIL_POW(sw,2));
    dZAA1 =  (5*TSIL_POW(e,2))/3.;
    dZZA1 =  (-4*TSIL_POW(e,2)*TSIL_POW(mw,2)*(TSIL_POW(cw,2) + TSIL_POW(sw,2)))/(cw*TSIL_POW(mz,2)*sw);
    dZAZ1 =  (TSIL_POW(e,2)*(TSIL_POW(cw,2)*(12*TSIL_POW(mw,2) + 19*TSIL_POW(mz,2)) + 12*TSIL_POW(mw,2)*TSIL_POW(sw,2) + TSIL_POW(mz,2)*(-8 + 9*TSIL_POW(sw,2))))/(3.*cw*TSIL_POW(mz,2)*sw);
    dZZZ1 =  (TSIL_POW(e,2)*(-17 + 39*TSIL_POW(cw,4) + 2*(16 + TSIL_POW(cw,2))*TSIL_POW(sw,2) - 17*TSIL_POW(sw,4)))/(12.*TSIL_POW(cw,2)*TSIL_POW(sw,2));
    dg2 =  4.0*TSIL_POW(e/sw,3);
    Ca =  -2.*sw;
    Cz1 =  1.;
    Cz2 =  -2.*cw;
    Cw1 =  -2.;
    Cw2 =  -2.;
    
    TSIL_REAL a = -1., b = 1., c = 1.;
    
    Aa = a*TSIL_A_ (ma2 , Q2);
    
    Ac = a*TSIL_A_ (MChi2 , Q2);
    
    Aw = a*TSIL_A_ (mw2 , Q2);
    
    Az = a*TSIL_A_ (mz2 , Q2);
    
    Bac = b*TSIL_B_ (ma2, MChi2, s, Q2);
    
    Bcw = b*TSIL_B_ (MChi2, mw2, s, Q2);
    
    Bcz = b*TSIL_B_ (MChi2, mz2, s, Q2);
    
    Bwc = Bcw;
    Bzc = Bcz;
    Bca = Bac;
    
    dBwc = c*TSIL_dBds_(mw2,MChi2,s,Q2);
    dBzc = c*TSIL_dBds_(mz2,MChi2,s,Q2);
    dBca = c*TSIL_dBds_(ma2,MChi2,s,Q2);
    
    dBac = c*TSIL_dBds_(ma2,MChi2,s,Q2);
    
  }
  
  double add_derivatives(Data &data)
  {
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*MChi*(-(((Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
                                                  (Ac*Power(MChi,2) - Aw*Power(MChi,2) - Power(MChi,4) - 4*dBwc*Power(MChi,6) + Bcw*Power(MChi,2)*Power(mw,2) -
                                                   2*dBwc*Power(MChi,4)*Power(mw,2)))/Power(MChi,6)) +
                                               ((2 + 8*dBwc*Power(MChi,2) + 8*Power(cw,2)*dBzc*Power(MChi,2) -
                                                 (2*(Ac - Aw - Bcw*Power(MChi,2) + 2*dBwc*Power(MChi,4) + Bcw*Power(mw,2) - dBwc*Power(MChi,2)*Power(mw,2)))/Power(MChi,2) +
                                                 (Ac - Aw - Bcw*(2*Power(MChi,2) - Power(mw,2)))/Power(MChi,2) -
                                                 (2*Power(cw,2)*(Ac - Az - Bcz*Power(MChi,2) + 2*dBzc*Power(MChi,4) + Bcz*Power(mz,2) - dBzc*Power(MChi,2)*Power(mz,2)))/
                                                 Power(MChi,2) + (Power(cw,2)*(Ac - Az - Bcz*(2*Power(MChi,2) - Power(mz,2))))/Power(MChi,2) +
                                                 8*dBac*Power(MChi,2)*Power(sw,2) - (2*(-Aa + Ac + Bca*Power(ma,2) - Bca*Power(MChi,2) - dBac*Power(ma,2)*Power(MChi,2) +
                                                                                        2*dBac*Power(MChi,4))*Power(sw,2))/Power(MChi,2) +
                                                 ((-Aa + Ac - Bca*(-Power(ma,2) + 2*Power(MChi,2)))*Power(sw,2))/Power(MChi,2))*
                                                (Aw + Az*Power(cw,2) + 2*Power(MChi,2) - 2*Bcw*Power(MChi,2) - 2*Bcz*Power(cw,2)*Power(MChi,2) - Bcw*Power(mw,2) -
                                                 Bcz*Power(cw,2)*Power(mz,2) + Aa*Power(sw,2) - Bca*Power(ma,2)*Power(sw,2) - 2*Bca*Power(MChi,2)*Power(sw,2) -
                                                 Ac*(1 + Power(cw,2) + Power(sw,2))))/(4.*Power(MChi,2))))/(64.*Power(Pi,4)*Power(sw,4));
    
    // These two expressions are symbolically identical yet will give slightly different (1e-20) numerical results
    /*
     TSIL_COMPLEXCPP result_2 = (Power(e,4)*(-4*(Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*(Ac - Aw - Power(MChi,2) - 4*dBwc*Power(MChi,4) + Bcw*Power(mw,2) - 2*dBwc*Power(MChi,2)*Power(mw,2)) +
     ((Aw*Power(cw,2) + Power(MChi,2) - 2*Bcz*Power(MChi,2) + Power(cw,2)*Power(MChi,2) - 2*Bcw*Power(cw,2)*Power(MChi,2) - Bcw*Power(cw,2)*Power(mw,2) - Bcz*Power(mz,2) + Aa*Power(cw,2)*Power(sw,2) -
     Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) - 2*Power(MChi,2)*Power(sw,2) + 4*Bcz*Power(MChi,2)*Power(sw,2) + Power(cw,2)*Power(MChi,2)*Power(sw,2) - 2*Bca*Power(cw,2)*Power(MChi,2)*Power(sw,2) +
     2*Bcz*Power(mz,2)*Power(sw,2) + Power(MChi,2)*Power(sw,4) - 2*Bcz*Power(MChi,2)*Power(sw,4) - Bcz*Power(mz,2)*Power(sw,4) + Az*Power(-1 + Power(sw,2),2) -
     Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2))))*(Aw*Power(cw,2) + Power(MChi,2) + Power(cw,2)*Power(MChi,2) + 4*Power(cw,2)*dBwc*Power(MChi,4) + 4*dBzc*Power(MChi,4) -
     Bcw*Power(cw,2)*Power(mw,2) + 2*Power(cw,2)*dBwc*Power(MChi,2)*Power(mw,2) - Bcz*Power(mz,2) + 2*dBzc*Power(MChi,2)*Power(mz,2) + Aa*Power(cw,2)*Power(sw,2) - Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) -
     2*Power(MChi,2)*Power(sw,2) + Power(cw,2)*Power(MChi,2)*Power(sw,2) + 2*Power(cw,2)*dBac*Power(ma,2)*Power(MChi,2)*Power(sw,2) + 4*Power(cw,2)*dBac*Power(MChi,4)*Power(sw,2) -
     8*dBzc*Power(MChi,4)*Power(sw,2) + 2*Bcz*Power(mz,2)*Power(sw,2) - 4*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,2) + Power(MChi,2)*Power(sw,4) + 4*dBzc*Power(MChi,4)*Power(sw,4) -
     Bcz*Power(mz,2)*Power(sw,4) + 2*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,4) + Az*Power(-1 + Power(sw,2),2) - Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2)))))/Power(cw,4)))/
     (256.*Power(MChi,3)*Power(Pi,4)*Power(sw,4));
     */
    
    
    return real(result);
    
  }
  
  double F11_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*(Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
                              (-Ac + Aw + Power(MChi,2) + 4*dBwc*Power(MChi,4) - Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2)))/
    (64.*Power(MChi,3)*Power(Pi,4)*Power(sw,4));
    
    return real(result);
    
  }
  
  
  
  double F12_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*(Aw*Power(cw,2) + Power(MChi,2) + Power(cw,2)*Power(MChi,2) + 4*Power(cw,2)*dBwc*Power(MChi,4) + 4*dBzc*Power(MChi,4) -
                                          Bcw*Power(cw,2)*Power(mw,2) + 2*Power(cw,2)*dBwc*Power(MChi,2)*Power(mw,2) - Bcz*Power(mz,2) + 2*dBzc*Power(MChi,2)*Power(mz,2) +
                                          Aa*Power(cw,2)*Power(sw,2) - Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) - 2*Power(MChi,2)*Power(sw,2) +
                                          Power(cw,2)*Power(MChi,2)*Power(sw,2) + 2*Power(cw,2)*dBac*Power(ma,2)*Power(MChi,2)*Power(sw,2) +
                                          4*Power(cw,2)*dBac*Power(MChi,4)*Power(sw,2) - 8*dBzc*Power(MChi,4)*Power(sw,2) + 2*Bcz*Power(mz,2)*Power(sw,2) -
                                          4*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,2) + Power(MChi,2)*Power(sw,4) + 4*dBzc*Power(MChi,4)*Power(sw,4) -
                                          Bcz*Power(mz,2)*Power(sw,4) + 2*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,4) + Az*Power(-1 + Power(sw,2),2) -
                                          Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2))))*
                              (-(Aw*Power(cw,2)) - Power(MChi,2) + 2*Bcz*Power(MChi,2) - Power(cw,2)*Power(MChi,2) + 2*Bcw*Power(cw,2)*Power(MChi,2) +
                               Bcw*Power(cw,2)*Power(mw,2) + Bcz*Power(mz,2) - Aa*Power(cw,2)*Power(sw,2) + Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) +
                               2*Power(MChi,2)*Power(sw,2) - 4*Bcz*Power(MChi,2)*Power(sw,2) - Power(cw,2)*Power(MChi,2)*Power(sw,2) +
                               2*Bca*Power(cw,2)*Power(MChi,2)*Power(sw,2) - 2*Bcz*Power(mz,2)*Power(sw,2) - Power(MChi,2)*Power(sw,4) +
                               2*Bcz*Power(MChi,2)*Power(sw,4) + Bcz*Power(mz,2)*Power(sw,4) - Az*Power(-1 + Power(sw,2),2) +
                               Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2)))))/(256.*Power(cw,4)*Power(MChi,3)*Power(Pi,4)*Power(sw,4));
    
    return real(result);
    
  }
  
  
  TSIL_COMPLEXCPP  gammagamma_Chi(Data data,double Q)
  {
    DoTSIL(data);
    
    p = Q;// TSIL_POW(data.Q,0.5);
    TSIL_REAL Q2 = pow(Q,2);
    
    TSIL_COMPLEXCPP AcMB = -i*TSIL_A_ (MChi2 ,  Q2);
    
    // evaluate as s = Q^2
    TSIL_COMPLEXCPP BccMB = i*TSIL_B_ (MChi2, MChi2, Q2, Q2);
    
    TSIL_COMPLEXCPP C0 = (4*Power(e,2)*(6*Power(MChi,2) - (-Power(p,2))))/9.;
    TSIL_COMPLEXCPP CAc =   Complex(0,2.6666666666666665)*Power(e,2) ;
    TSIL_COMPLEXCPP CBcc =   Complex(0,-1.3333333333333333)*Power(e,2)*(2*Power(MChi,2) + Power(p,2)) ;
    
    TSIL_COMPLEXCPP result = + C0  + AcMB * CAc + BccMB * CBcc;
    
    return -result/(16.0L*TSIL_POW(PI,2));
  }
  
  
  
  
}



// determine MSbar parameters
void set_SM_parameters_2loop(Data &data)
{
  Self_energy se;
  
  
  
  // matching (?) of SM to Wino model
  data.alpha = data.alpha*( 1.0-real(extra_TSIL_interface::gammagamma_Chi(data,TSIL_POW(data.Q,0.5))) /data.Q );
  
  // determine MS bar input parameters at Q
  // only recompute the 1-loop functions for efficiency
  data.do_tsil_all = false;
  
  data.P = data.mw;
  se.run_tsil(data);
  data.mw = pow( pow(data.mw,2) - real(data.SE_1["V3"]),0.5 );
  
  data.P = data.mz;
  se.run_tsil(data);
  data.mz = pow( pow(data.mz,2) - real(data.SE_1["V2"]) ,0.5);
  
  cout << "mw, mz = " << data.mw << " , " << data.mz << endl;
  
  // need to iterate to determine MS bar mass for MChi to match equation (9) of Ibe et al.
  //data.MChi = iterative_ms_bar_mass(data, "F11_g1");
  
  data.do_tsil_all = true;
}



// determine MSbar parameters
void set_SM_parameters_1loop(Data &data)
{
  Self_energy se;
  
  
  // RGE evolution
  double A = (19./(10.*Pi)); // MSSM
  //double A = (17./(30.*Pi)); // Wino model
  //double A = (7./(30.*Pi)); // SM
  double alpha_mz = data.alpha;
  double mu0 = data.mz;
  double mu = pow(data.Q,0.5);
  
  data.alpha = pow(   1.0/alpha_mz -  A * log( mu/mu0) , -1);
  
}



double iterative_mass(Data data, string particle,int loop_order)
{
  long double M_tree = data.MChi;
  long double M_pole = data.MChi;
  data.P = M_tree;
  long double diff = 1;
  long double precision = 1e-5;
  int iteration =0;
  int loop_order_temp = loop_order;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    if (iteration==0)
    {
     loop_order = 1;
    }
    
    if (loop_order == 1)
    {
      M_pole = M_tree + (data.SE_1[particle]);
    }
    
    if (loop_order == 2)
    {
      M_pole = M_tree + (data.SE_1[particle]+data.SE_2[particle]);
    }
    
    diff = abs(M_pole - data.P);
    
    //diff = abs(M_pole/data.P - 1.0);
    
    data.P = M_pole;
    
    
    if (iteration == 0)
    {
      loop_order = loop_order_temp;
    }
    
    iteration++;
    //cout<< "\r" << "M_pole - p = " << diff << " GeV";
    //std::cout << std::flush;
  } while (diff > precision  && iteration < 50000);
  
  //cout<< "\r" << "M_pole - p = " << diff << " GeV";
  //cout << "\n";
  
  
  if (iteration == 50000)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  cout << "number of iterations = " << iteration << endl;
  
  return M_pole;
}



void plot_M(Data data)
{
  
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 100;
  double n = 0;
  double M = 0;
  int status = 0;
  
  double max_M = 1e6; // (GeV)
  double min_M = 10; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  double alpha_in = data.alpha;
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    M = pow(10,n);
    data.MChi=M;
    data.P = M;
    data.alpha = alpha_in;
    data.mz = mz_in;
    data.mw = mw_in;
    
    
    data.Q = pow(2*M,2);
    
    // compute explicit mass splitting
    set_SM_parameters_1loop(data);
    data.P = M;
    Self_energy se;
    se.run_tsil(data);
    
    //double delta_m_explicit_2loop = (data.SE_1["F12_g1"]+data.SE_2["F12_g1"]) - (data.SE_1["F11_g1"]+data.SE_2["F11_g1"]) - extra_TSIL_interface::add_derivatives(data);
    
    double delta_m_explicit_1loop = (data.SE_1["F12_g1"]) - (data.SE_1["F11_g1"]);
    
    // compute iterative mass splitting
    
    double it_1loop_F11 = iterative_mass(data,"F11_g1",1);
    double it_1loop_F12 = iterative_mass(data,"F12_g1",1);
    
    
    //double delta_m_iterative_2loop = iterative_mass(data,"F12_g1",2) - iterative_mass(data,"F11_g1",2);
    
    double delta_m_iterative_1loop = it_1loop_F12 - it_1loop_F11;
    
    
    myfile << M << " " << delta_m_explicit_1loop <<  " " << delta_m_iterative_1loop << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  cout << "example mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_itMSSM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  myfile.close();
}



void print_masses(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  
  double F11_explicit_2loop = data.MChi + data.SE_1["F11_g1"] + data.SE_2["F11_g1"] + extra_TSIL_interface::F11_der(data);
  
  double F12_explicit_2loop = data.MChi + data.SE_1["F12_g1"] + data.SE_2["F12_g1"] + extra_TSIL_interface::F12_der(data);
  
  double F12_explicit_1loop = data.MChi + data.SE_1["F12_g1"];
  
  double F11_explicit_1loop = data.MChi + data.SE_1["F11_g1"];
  
  
  cout << "Self energies" << endl;
  
  cout << "One-loop self energy (charged) " <<data.SE_1["F12_g1"] << endl;
  cout << "One-loop self energy (neutral) " <<data.SE_1["F11_g1"] << endl;
  
  
  double F11_iterative_2loop = 0;//iterative_mass(data, "F11_g1" ,2);
  
  double F12_iterative_2loop = 0;//iterative_mass(data, "F12_g1" ,2);
  
  data.do_tsil_all = false;
  
  double F11_iterative_1loop = iterative_mass(data, "F11_g1" ,1);
  
  double F12_iterative_1loop = iterative_mass(data, "F12_g1" ,1);
  
  cout << "2-loop pole masses" << endl;
  cout << "iterative pole mass (neutral) = " << F11_iterative_2loop << endl;
  cout << "explicit pole mass (neutral) = " << F11_explicit_2loop << endl;
  
  cout << "iterative pole mass (charged)  = " << F12_iterative_2loop << endl;
  cout << "explicit pole mass (charged) = " << F12_explicit_2loop << endl;
  
  cout << "Mass-splitting" << endl;
  
  cout << "1-loop iterative mass splitting = " << F12_iterative_1loop - F11_iterative_1loop << endl;
  
  cout << "2-loop iterative mass splitting = " << F12_iterative_2loop - F11_iterative_2loop << endl;
  
  cout << "1-loop explicit mass splitting = " << F12_explicit_1loop - F11_explicit_1loop << endl;
  
  cout << "2-loop explicit mass splitting = " << F12_explicit_2loop - F11_explicit_2loop << endl;
  
}








int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  Data data(options);
  
  

  //print_masses(data);
  
  plot_M(data);
  
  
  
  
  return 0;
}
