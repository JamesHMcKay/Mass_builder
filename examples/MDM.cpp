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
  
  long double  MChi,  MChi2 ,  ma,  ma2 ,  mw,  mw2 ,  mz,  mz2 ,  mh,  mh2 ,  null,  null2 ,  mf,  mf2 ,  mt, mt2 ;
  
  long double Pi;
  long double CTW, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, STW, cw, cw2, d5Z, d5m, d6Z, d6m, d7Z, d7m, dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, g, g1, g2, sw, sw2, v ;
  
  
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca,Bac,Bcz,Bcw;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,dBac,M;
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az;
  
  
  void DoTSIL(Data data)
  {
    
    
    cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   STW = data.STW,   CTW = data.CTW,   g2 = data.g2,   g1 = data.g1,   g = data.g,   v = data.v ;
    
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
    g2 =  2.*mw/v;
    g1 =  g2*sw/cw;
    g =  g2;
    
    
    
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
  
  
  double F5_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(g2,4)*(5*Aw + Az*Power(CTW,2) + 5*Power(MChi,2) + Power(CTW,2)*Power(MChi,2) + 20*dBwc*Power(MChi,4) + 
       4*Power(CTW,2)*dBzc*Power(MChi,4) - 5*Bcw*Power(mw,2) + 10*dBwc*Power(MChi,2)*Power(mw,2) - Bcz*Power(CTW,2)*Power(mz,2) + 
       2*Power(CTW,2)*dBzc*Power(MChi,2)*Power(mz,2) + Aa*Power(STW,2) - Bca*Power(ma,2)*Power(STW,2) + Power(MChi,2)*Power(STW,2) + 
       2*dBac*Power(ma,2)*Power(MChi,2)*Power(STW,2) + 4*dBac*Power(MChi,4)*Power(STW,2) - Ac*(5 + Power(CTW,2) + Power(STW,2)))*
     (-5*Aw - Az*Power(CTW,2) - 5*Power(MChi,2) + 10*Bcw*Power(MChi,2) - Power(CTW,2)*Power(MChi,2) + 2*Bcz*Power(CTW,2)*Power(MChi,2) + 
       5*Bcw*Power(mw,2) + Bcz*Power(CTW,2)*Power(mz,2) - Aa*Power(STW,2) + Bca*Power(ma,2)*Power(STW,2) - Power(MChi,2)*Power(STW,2) + 
       2*Bca*Power(MChi,2)*Power(STW,2) + Ac*(5 + Power(CTW,2) + Power(STW,2))))/(256.*Power(MChi,3)*Power(Pi,4));
    
    return real(result);
    
  }
  
  
  
  double F6_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(g2,4)*(Ac - Aw + 2*Ac*Power(CTW,2) - 2*Az*Power(CTW,2) - Power(MChi,2) + 2*Bcw*Power(MChi,2) - 2*Power(CTW,2)*Power(MChi,2) + 
       4*Bcz*Power(CTW,2)*Power(MChi,2) + Bcw*Power(mw,2) + 2*Bcz*Power(CTW,2)*Power(mz,2) - 2*Aa*Power(STW,2) + 2*Ac*Power(STW,2) + 
       2*Bca*Power(ma,2)*Power(STW,2) - 2*Power(MChi,2)*Power(STW,2) + 4*Bca*Power(MChi,2)*Power(STW,2))*
     (Aw + 2*Az*Power(CTW,2) + Power(MChi,2) + 2*Power(CTW,2)*Power(MChi,2) + 4*dBwc*Power(MChi,4) + 8*Power(CTW,2)*dBzc*Power(MChi,4) - 
       Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2) - 2*Bcz*Power(CTW,2)*Power(mz,2) + 
       4*Power(CTW,2)*dBzc*Power(MChi,2)*Power(mz,2) + 2*Aa*Power(STW,2) - 2*Bca*Power(ma,2)*Power(STW,2) + 
       2*Power(MChi,2)*Power(STW,2) + 4*dBac*Power(ma,2)*Power(MChi,2)*Power(STW,2) + 8*dBac*Power(MChi,4)*Power(STW,2) - 
       Ac*(1 + 2*Power(CTW,2) + 2*Power(STW,2))))/(64.*Power(MChi,3)*Power(Pi,4));
    return real(result);
    
  }
  
  double F7_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (9*Power(g2,4)*(Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
     (-Ac + Aw + Power(MChi,2) + 4*dBwc*Power(MChi,4) - Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2)))/
   (64.*Power(MChi,3)*Power(Pi,4));
    
    return real(result);
    
  }
  
}




double iterative_mass(Data data, string particle,int loop_order)
{
  double M_tree = data.MChi;
  double M_pole = data.MChi;
  data.P = M_tree;
  double diff = 1;
  double precision = 1e-4;
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
    
    data.P = M_pole;
    
    if (iteration == 0)
    {
      loop_order = loop_order_temp;
    }
    
    iteration++;
    cout<< "\r" << "M_pole - p = " << diff << " GeV";
    std::cout << std::flush;
  } while (diff > precision  && iteration < 500);
  
  cout<< "\r" << "M_pole - p = " << diff << " GeV";
  cout << "\n";
  
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  return M_pole;
}








void plot_M(Data data)
{
  
  
  ofstream myfile;
  myfile.open ("models/MDM/output/mass_splittings.txt");
  int pts = 10;
  double n = 0;
  double M = 0;
  int status = 0;
  
  double max_M = 10000; // (GeV)
  double min_M = 10; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    
    
    M = pow(10,n);
    data.MChi=M;
    data.P = M;
    
    //data.Q = M;
    
    // compute explicit mass splitting
    
    Self_energy se;
    se.run_tsil(data);
    
  
    double delta_m1_explicit_1loop = (data.SE_1["F5"]) - (data.SE_1["F7"]);
    double delta_m2_explicit_1loop = (data.SE_1["F6"]) - (data.SE_1["F7"]);
    
    double delta_m1_explicit_2loop = (data.SE_1["F5"]+data.SE_2["F5"]+extra_TSIL_interface::F5_der(data)) - (data.SE_1["F7"]+data.SE_2["F7"]+extra_TSIL_interface::F7_der(data));
    double delta_m2_explicit_2loop = (data.SE_1["F6"]+data.SE_2["F6"]+extra_TSIL_interface::F6_der(data)) - (data.SE_1["F7"]+data.SE_2["F7"]+extra_TSIL_interface::F7_der(data));
    
    // compute iterative mass splitting
    
    
    
    //double delta_m1_iterative_1loop = iterative_mass(data,"F5",1) - iterative_mass(data,"F7",1);
    
    //double delta_m2_iterative_1loop = iterative_mass(data,"F6",1) - iterative_mass(data,"F7",1);
    
    
    myfile << M << " " << delta_m1_explicit_2loop  <<  " " << delta_m2_explicit_2loop << " " << delta_m1_explicit_1loop <<  " " << delta_m2_explicit_1loop << endl;
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
  
  return 0;
}
