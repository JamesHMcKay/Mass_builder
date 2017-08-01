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
  long double  MChi,  MChi2 ,  ma,  ma2 ,  mw,  mw2 ,  mz,  mz2 ,  mh,  mh2 ,  null, null2 ;
  
  double Pi;
  long double CTW, Ca, Cw1, Cw2, Cz1, Cz2, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, MassBuilderJEpsilon, STW, cw, cw2, d1Z, d1m, d2Z, d2m, dMWsq1, dMZsq1, dZAA1, dZAZ1, dZW1, dZZA1, dZZZ1, dg2, e, g, g1, g2, sw, sw2, v ;
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca,Bac,Bcz,Bcw;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,dBac,M;
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az;
  
  
  void DoTSIL(Data data)
  {
    cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   STW = data.STW,   CTW = data.CTW,   g1 = data.g1,   g2 = data.g2,   g = data.g,   d1Z = data.d1Z,   d1m = data.d1m,   d2Z = data.d2Z,   d2m = data.d2m,   dZW1 = data.dZW1,   dMWsq1 = data.dMWsq1,   dMZsq1 = data.dMZsq1,   dZAA1 = data.dZAA1,   dZZA1 = data.dZZA1,   dZAZ1 = data.dZAZ1,   dZZZ1 = data.dZZZ1,   dg2 = data.dg2,   e = data.e,   Ca = data.Ca,   Cz1 = data.Cz1,   Cz2 = data.Cz2,   Cw1 = data.Cw1,   Cw2 = data.Cw2,   v = data.v,   MassBuilderJEpsilon = data.MassBuilderJEpsilon ;
    
    TSIL_REAL Q2 = data.Q;
    TSIL_REAL s = pow(data.P,2);
    MassBuilderCTM1 = 0;
    MassBuilderCTZ1 = 0;
    MassBuilderCTM2 = 0;
    MassBuilderCTZ1 = 0;
    MChi = data.MChi, MChi2 = TSIL_POW(data.MChi, 2) ,   ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,   mw = data.mw, mw2 = TSIL_POW(data.mw, 2) ,   mz = data.mz, mz2 = TSIL_POW(data.mz, 2) ,   mh = data.mh, mh2 = TSIL_POW(data.mh, 2) ,   null = data.null , null2 = TSIL_POW(data.null, 2) ;
    
    dcomp ii=-1;ii=sqrt(ii);i=ii;
    Pi=PI;
    cw =  mw/mz;
    cw2 =  TSIL_POW(cw,2);
    sw =  TSIL_POW(1.-cw2,0.5);
    sw2 =  TSIL_POW(sw,2);
    STW =  sw;
    CTW =  cw;
    g1 =  2.*mw/v;
    g2 =  g1*sw/cw;
    g =  g2;
    e =  g1*g2/TSIL_POW(g1*g1+g2*g2,0.5);
    d1Z =  (-2*TSIL_POW(e,2)*(TSIL_POW(-1 + TSIL_POW(sw,2),2) + TSIL_POW(cw,2)*(1 + TSIL_POW(sw,2))))/(TSIL_POW(cw,2)*TSIL_POW(sw,2));
    d1m =  (8*TSIL_POW(e,2)*MChi*(TSIL_POW(-1 + TSIL_POW(sw,2),2) + TSIL_POW(cw,2)*(1 + TSIL_POW(sw,2))))/(TSIL_POW(cw,2)*TSIL_POW(sw,2));
    d2Z =  (-4*TSIL_POW(e,2))/TSIL_POW(sw,2);
    d2m =  (16*TSIL_POW(e,2)*MChi)/TSIL_POW(sw,2);
    dZW1 =  -(TSIL_POW(e,2)*(-9 + 20*TSIL_POW(cw,2) + 20*TSIL_POW(sw,2)))/(6.*TSIL_POW(sw,2));
    dMWsq1 =  (TSIL_POW(e,2)*(TSIL_POW(cw,4)*(44*TSIL_POW(mw,2) + 6*TSIL_POW(mz,2)) - 6*TSIL_POW(mw,2)*TSIL_POW(sw,4) + TSIL_POW(cw,2)*(6*TSIL_POW(ma,2)*TSIL_POW(sw,2) + TSIL_POW(mw,2)*(-33 + 38*TSIL_POW(sw,2)))))/(6.*TSIL_POW(cw,2)*TSIL_POW(sw,2));
    dMZsq1 =  (TSIL_POW(e,2)*(TSIL_POW(cw,2)*TSIL_POW(mz,2)*(-17 + 39*TSIL_POW(cw,4) + 2*(16 + TSIL_POW(cw,2))*TSIL_POW(sw,2) - 17*TSIL_POW(sw,4)) +  12*TSIL_POW(mw,2)*(-1 + 2*TSIL_POW(cw,6) - 2*TSIL_POW(cw,2)*TSIL_POW(sw,4))))/(12.*TSIL_POW(cw,4)*TSIL_POW(sw,2));
    dZAA1 =  (-5*TSIL_POW(e,2))/3.;
    dZZA1 =  (4*TSIL_POW(e,2)*TSIL_POW(mw,2)*(TSIL_POW(cw,2) + TSIL_POW(sw,2)))/(cw*TSIL_POW(mz,2)*sw);
    dZAZ1 =  -(TSIL_POW(e,2)*(TSIL_POW(cw,2)*(12*TSIL_POW(mw,2) + 19*TSIL_POW(mz,2)) + 12*TSIL_POW(mw,2)*TSIL_POW(sw,2) + TSIL_POW(mz,2)*(-8 + 9*TSIL_POW(sw,2))))/(3.*cw*TSIL_POW(mz,2)*sw);
    dZZZ1 =  (TSIL_POW(e,2)*(17 - 39*TSIL_POW(cw,4) - 2*(16 + TSIL_POW(cw,2))*TSIL_POW(sw,2) + 17*TSIL_POW(sw,4)))/(12.*TSIL_POW(cw,2)*TSIL_POW(sw,2));
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
    
    
    return real(result);
  
  }
  
}


double iterative_mass_F5(Data data)
{
  
  double M = data.MChi;
  double Mp = data.MChi;
  
  double MFn=M,M_tree=M,new_MFn,old_MFn=M,p;
  
  double diff = 1;
  double precision = 0.0001;
  int iteration =0;
  
  do{
    p=old_MFn;
    
    data.P = p;
    Self_energy se;
    se.run_tsil(data);
    
    double M_1loop=M_tree + (data.SE_1["F12_g1"]+data.SE_2["F12_g1"]);
    MFn=M_1loop;
    new_MFn=MFn;
    diff=abs(new_MFn-old_MFn);
    old_MFn=new_MFn;
    iteration++;
    //cout << "diff = " << diff << endl;
    
  } while (diff > precision  && iteration < 500);
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  Mp = new_MFn;
  
  //cout << "----- done ----- " << endl;
  return Mp;
}

double iterative_mass_F6(Data data)
{
  
  double M = data.MChi;
  double Mp = data.MChi;
  
  double MFn=M,M_tree=M,new_MFn,old_MFn=M,p;
  
  double diff = 1;
  double precision = 0.0001;
  int iteration =0;
  
  do{
    p=old_MFn;
    
    data.P = p;
    Self_energy se;
    se.run_tsil(data);
    
    double M_1loop=M_tree +(data.SE_1["F11_g1"]+data.SE_2["F11_g1"]);
    MFn=M_1loop;
    new_MFn=MFn;
    diff=abs(new_MFn-old_MFn);
    old_MFn=new_MFn;
    iteration++;
    //cout << "diff = " << diff << endl;
    
  } while (diff > precision && iteration < 500);
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  //cout << "----- done ----- " << endl;
  
  Mp = new_MFn;
  return Mp;
}


int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  
  Data data(options);
  
  //cout << " derivative terms = " << der << endl;
  //cout << " two-loop mass splitting = " << pole_mass_F6(data) - pole_mass_F5(data)+der << endl;
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 10;
  double n = 0;
  double M = 0;
  int status = 0;
  for (int i = 0; i < pts ; i++)
  {
    n=(float(i)/float(pts))*5;
    
    // M= pow(10,n);
    // data.MChi=M;
    // data.P = M;
    //double delta_m_it=iterative_mass_F5(data) - iterative_mass_F6(data);
    M = pow(10,n);
    data.MChi=M;
    data.P = M;
    
    Self_energy se;
    se.run_tsil(data);
    
    double delta_m_2 = data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - extra_TSIL_interface::add_derivatives(data);
    
    double delta_m_1 = data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
    
    myfile << M << " " << delta_m_1 <<  " " << delta_m_1 + delta_m_2 << endl;
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
  
  
  return 0;
}
