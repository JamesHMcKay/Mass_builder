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

double Pi=PI;

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



double iterative_ms_bar_mass(Data data, string particle)
{
  double M_msbar = data.MChi;
  double M_pole = data.MChi;
  data.P = M_msbar;
  double diff = 1;
  double precision = 1e-6;
  int iteration =0;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    
    M_msbar = M_pole - (data.SE_1[particle]);
    
    diff = abs(M_msbar - data.MChi);
    
    data.MChi = M_msbar;
    
    iteration++;
    // cout<< "\r" << "M_msbar - MChi = " << diff << " GeV";
    // std::cout << std::flush;
  } while (diff > precision  && iteration < 500);
  
  // cout<< "\r" << "M_msbar - MChi = " << diff << " GeV";
  // cout << "\n";
  
  
  if (iteration == 500)
  {
    cout << "ms bar mass did not converge" << endl;
  }
  
  return M_msbar;
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
  
  // need to iterate to determine MS bar mass for MChi to match equation (9) of Ibe et al.
  //data.MChi = iterative_ms_bar_mass(data, "F11_g1");
  
  data.do_tsil_all = true;
}


void plot_Q(Data data)
{
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings_Q.txt");
  int pts = 50;
  double M = 1000;
  int status = 0;
  
  double max_Q = 400; // (GeV)
  double min_Q = 50; // (GeV)
  data.MChi=M;
  data.P = M;
  
  double alpha_in = data.alpha;
  double MChi_in = data.MChi;
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  std::vector<double> mt={ 176.424,
    174.888,
    173.566,
    172.409,
    171.382,
    170.462,
    169.629,
    168.868,
    168.17,
    167.525,
    166.926,
    166.367,
    165.844,
    165.353,
    164.89,
    164.452,
    164.037,
    163.643,
    163.268,
    162.91,
    162.568,
    162.24,
    161.926,
    161.625,
    161.335,
    161.055,
    160.786,
    160.526,
    160.275,
    160.033,
    159.798,
    159.57,
    159.349,
    159.135,
    158.927,
    158.725,
    158.528,
    158.337,
    158.15,
    157.969,
    157.792,
    157.619,
    157.451,
    157.286,
    157.125,
    156.968,
    156.814,
    156.664,
    156.517,
    156.372,
    156.231};
  
  
  for (int i = 0; i < pts+1 ; i++)
  {
    
    data.mt = mt[i];
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    
    
    set_SM_parameters_1loop(data);
    data.P = data.MChi;
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    data.do_tsil_all = false;
    Self_energy se1;
    se1.run_tsil(data);
    data.do_tsil_all = true;
    
    double delta_m_1 = data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
    
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.P = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    set_SM_parameters_2loop(data);
    data.P = data.MChi;
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    Self_energy se2;
    se2.run_tsil(data);
    
    
    double delta_m_2 = (data.SE_1["F12_g1"] - data.SE_1["F11_g1"]) + ( data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - extra_TSIL_interface::add_derivatives(data) );
    
    myfile << pow(data.Q,0.5) << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  cout << "mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_MSSM_Q.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  myfile.close();
  
}


void plot_M(Data data)
{
  
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 10;
  double n = 0;
  double M = 0;
  int status = 0;
  
  double max_M = 4000; // (GeV)
  double min_M = 90; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
      
  double alpha_in = data.alpha;
  
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    
    M = pow(10,n);
    double MChi_in = M;
    data.MChi=M;
    data.P = M;

    
    
    data.mt = 163.3;
    
    data.Q = pow(163.3,2);
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    
    
    set_SM_parameters_1loop(data);
    data.P = data.MChi;
    data.Q = pow(163.3,2);
    
    data.do_tsil_all = false;
    Self_energy se1;
    se1.run_tsil(data);
    data.do_tsil_all = true;
    
    double delta_m_1 = data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
    
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.P = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    set_SM_parameters_2loop(data);
    data.P = data.MChi;
    data.Q = pow(163.3,2);
    
    Self_energy se2;
    se2.run_tsil(data);
    
    
    double delta_m_2 = (data.SE_1["F12_g1"] - data.SE_1["F11_g1"]) + ( data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - extra_TSIL_interface::add_derivatives(data) );
    
    myfile << M << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  cout << "example mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_MSSM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  myfile.close();
}


int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  
  Data data(options);
  
  plot_M(data);
  /*
  
  data.MChi = 1000;
  
  double alpha_in = data.alpha;
  double MChi_in = data.MChi;
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  double Qin = data.Q;
  
  data.alpha = alpha_in;
  data.MChi = MChi_in;
  data.P = MChi_in;
  data.mw = mw_in;
  data.mz = mz_in;
  set_SM_parameters_2loop(data);
  data.P = data.MChi;
  data.Q = Qin;
  
  Self_energy se2;
  se2.run_tsil(data);
  
  
  cout << "two-loop mass splitting = " << (data.SE_1["F12_g1"] - data.SE_1["F11_g1"]) + ( data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - extra_TSIL_interface::add_derivatives(data) ) << endl;
  */
  
  return 0;
}
