/*
 Mass Builder 
 
 James McKay
 Sep 2016
 
 -- EW_triplet.cpp --
 
 compute full two-loop self energy including derivative of 1-loop functions
 
 requires an input list flag at runtime: ./EW_triplet -i models/EW_triplet/input.txt
 */

#include "data.hpp"
#include "self_energy.hpp"
#include "EW_triplet.hpp"
#include "cmake_variables.hpp"



using namespace std;


// extra TSIL interface for manually entered derivatives
// of one-loop self energies

namespace extra_TSIL_interface
{
  #include TSIL_PATH
  using namespace std;
  TSIL_DATA bar;
  
  #ifndef PI
  #define PI 4.0L*atan(1.0L)
  #endif
  long double strtold(const char *, char **);
  // define subroutines here
  TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}
  TSIL_REAL Sin(TSIL_REAL a){return sin(a);}
  TSIL_REAL Cos(TSIL_REAL a){return cos(a);}
  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}
  TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az, Bac ;
  TSIL_COMPLEXCPP Bcw, Bcz ;
  TSIL_COMPLEXCPP  i=Power(-1,0.5);
  TSIL_REAL ma, mc, mw, mz ;
  TSIL_REAL ma2,mc2,mw2,mz2;
  
  TSIL_REAL Pi,g,g1, tW, sw2, cw2, sw, cw , C, p, v;
  
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,M;
  
  
  TSIL_REAL dBds(double m, double ma, double p)
  {
    m = p;
    return (1.0/m)+ (1.0/(2.0*m))*log( ma/m );
  }
  
  void DoTSIL_2(TSIL_REAL s,TSIL_REAL Q2)
  {
    Aa = -i*TSIL_A_ (ma2 , Q2);
    
    Ac = -i*TSIL_A_ (mc2 , Q2);
    
    Aw = -i*TSIL_A_ (mw2 , Q2);
    
    Az = -i*TSIL_A_ (mz2 , Q2);
    
    Bac = i*TSIL_B_ (ma2, mc2, s, Q2);
    
    Bcw = i*TSIL_B_ (mc2, mw2, s, Q2);
    
    Bcz = i*TSIL_B_ (mc2, mz2, s, Q2);
    
    Bwc = Bcw;
    Bzc = Bcz;
    Bca = Bac;
    
    dBwc = i*TSIL_dBds_(mw2,mc2,s,Q2);
    dBzc = i*TSIL_dBds_(mz2,mc2,s,Q2);
    dBca = i*TSIL_dBds_(ma2,mc2,s,Q2);
  }
  
  int init(Data data)
  {
  
    ;
  
    cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   g = data.g,   v = data.v;
  
  
    mw= data.mw, mz = data.mz ,ma = data.ma, mc = data.MChi;
    
    cw =  mw/mz;
    cw2 =  TSIL_POW(cw,2);
    sw =  TSIL_POW(1.-cw2,0.5);
    
    sw2 =  TSIL_POW(sw,2);
    g1 =  2.*mw/v;
    g =  g1*sw/cw;
    
    mw2 = pow(mw,2);
    mz2 = pow(mz,2);
    ma2 = pow(ma,2);
    mc2 = pow(mc,2);
    tW = acos(mw/mz);
    sw2 = data.sw2, cw2 = data.cw2;
    sw = data.sw, cw = data.cw;
    C = TSIL_POW(g,2)/(16.0L*PI*PI);
    i=Power(-1,0.5);
    Pi=PI;
    DoTSIL_2( TSIL_POW(data.P,2),data.Q);
    return 0;
  }
  
  TSIL_COMPLEXCPP  SigmaK_0(TSIL_REAL p, TSIL_REAL M)
  {
    return (-2.0*C/TSIL_POW(p,2))*( (TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mw,2))*Bwc-Ac+Aw);
  }
  TSIL_COMPLEXCPP  SigmaK_1(TSIL_REAL p, TSIL_REAL M)
  {
    return (C/TSIL_POW(p,2)) * (cw2*( -(TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mz,2))*Bzc-Az )-sw2*(TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(ma,2))*Bca - (TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mw,2))*Bwc-sw2*Aa+2.0L*Ac-Aw);
  }
  
  TSIL_COMPLEXCPP  SigmaM_1(TSIL_REAL M)
  {
    return 4.0*M*C*(Bwc+cw2*Bzc+sw2*Bca);
  }
  
  TSIL_COMPLEXCPP  SigmaM_0(TSIL_REAL M)
  {
    return 8.0*M*C*Bwc;
  }
  
  // the following are derivatives of the one loop self energy functions
  TSIL_COMPLEXCPP  d_SigmaM_1(TSIL_REAL M)
  {
    return 4.0*M*C*(dBwc+cw2*dBzc+sw2*dBca);
  }
  
  TSIL_COMPLEXCPP  d_SigmaK_1(TSIL_REAL p, TSIL_REAL M)
  {
    return -SigmaK_1(p,M)/(TSIL_POW(p,2)) + (C/TSIL_POW(p,2)) *(-cw2*Bzc-sw2*Bca-Bwc+cw2*( -(TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mz,2))*dBzc )-sw2*(TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(ma,2))*dBca - (TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mw,2))*dBwc );
    
  }
  
  TSIL_COMPLEXCPP  d_SigmaM_0(TSIL_REAL M)
  {
    return 8.0*M*C*dBwc;
  }
  
  TSIL_COMPLEXCPP  d_SigmaK_0(TSIL_REAL p, TSIL_REAL M)
  {
    return -SigmaK_0(p,M)/(TSIL_POW(p,2))-(2.0*C/TSIL_POW(p,2))*(Bwc+(TSIL_POW(p,2)+TSIL_POW(M,2)-TSIL_POW(mw,2))*dBwc);
  }
  
  
  double One_loop_derivatives::add_derivatives(Data &data)
  {
    TSIL_REAL p = data.P, M = data.MChi;
    init(data);
    TSIL_COMPLEXCPP SE_0 = (SigmaM_0(M) + M*SigmaK_0(M,p) )*( SigmaK_0(M,p) + 2.0L*TSIL_POW(M,2)*d_SigmaK_0(M,p)+2.0L*M*d_SigmaM_0(M));
    TSIL_COMPLEXCPP SE_1 =(SigmaM_1(M) + M*SigmaK_1(M,p) )*( SigmaK_1(M,p) + 2.0L*TSIL_POW(M,2)*d_SigmaK_1(M,p)+2.0L*M*d_SigmaM_1(M));
    
    return real(SE_1-SE_0);
    
  }
  
  
  double One_loop_derivatives::add_derivatives_2(Data &data)
  {
    
  

    init(data);
    
    cout << "256.*Power(cw,4)*Power(Pi,4)*Power(sw,4));" << 256.*Power(cw,4)*Power(Pi,4)*Power(sw,4) << endl;
    
    TSIL_COMPLEXCPP result = (-3*Power(e,4)*MChi*Power(-1 + Power(sw,2),2)*(-3*Bwc*Bwc + Complex(0,8)*(Bwc - Bzc) + 2*Bwc*Bzc +Bzc*Bzc - Complex(0,12)*Power(MChi,2)*dBwc +
       18*Power(MChi,2)*Bwc*dBwc - 6*Power(MChi,2)*Bzc*dBwc + Complex(0,12)*Power(MChi,2)*dBzc - 
       6*Power(MChi,2)*Bwc*dBzc - 6*Power(MChi,2)*Bzc*dBzc + 
       Power(sw,4)*(Bac - Bzc)*(Bac - Bzc + 6*Power(MChi,2)*(-dBca + dBzc)) + 
       2*Power(sw,2)*(-Bzc*Bzc + 3*Power(MChi,2)*(Complex(0,-2) + Bwc)*(-dBca + dBzc) + 
          Bac*(Complex(0,-4) + Bwc + Bzc - 3*Power(MChi,2)*(dBwc + dBzc)) + 
          Bzc*(Complex(0,4) - Bwc + 3*Power(MChi,2)*(-dBca + dBwc + 2*dBzc)))))/(256.*Power(cw,4)*Power(Pi,4)*Power(sw,4));
  
    return real(result);
  
  }
  
  
  
}



Self_energy self_energy;

double pole_mass_F1(Data data)
{
  double Mp = data.MChi + (data.SE_1["F1"]+data.SE_2["F1"]);
  return Mp;
}

double pole_mass_F2(Data data)
{
  double Mp = data.MChi + (data.SE_1["F2"]+data.SE_2["F2"]);
  return Mp;
}

using namespace extra_TSIL_interface;


int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy self_energy;
  Data data(options);

  self_energy.run_tsil(data);
  
  cout << "one-loop mass splitting = " << data.SE_1["F11_g1"] - data.SE_1["F12_g1"] << endl;

  cout << "two-loop mass splitting = " << data.SE_2["F11_g1"] - data.SE_2["F12_g1"] << endl;
  
  
  
  One_loop_derivatives one_loop_derivatives(data);
  

  cout << "one-loop derivative terms = " << one_loop_derivatives.add_derivatives(data) << endl;
  cout << "one-loop derivative terms (alternative) = " << one_loop_derivatives.add_derivatives_2(data) << endl;


  cout << "two-loop mass splitting including derivative terms = " << data.SE_2["F11_g1"] - data.SE_2["F12_g1"]+one_loop_derivatives.add_derivatives(data) << endl;
  
  return 0;
}
