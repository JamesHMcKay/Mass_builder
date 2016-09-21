#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>

#include "supplements.hpp"


namespace supplementary_code
{
  /*TSIL_INCLUDE_PATH */#include "/Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h"  // Required TSIL header file
  
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
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az, Bac ;
  TSIL_COMPLEXCPP Bcw, Bcz ;
  TSIL_COMPLEXCPP  i=Power(-1,0.5);
  TSIL_REAL ma, mc, mw, mz ;
  
  TSIL_REAL Pi,g2, tW, sw2, cw2, sw, cw , C,p,S2TW;
  
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca;
  
  void DoTSIL_2(TSIL_REAL s,TSIL_REAL Q2)
  {
    cout << "mass of photon is = " << ma << endl;
    Aa = -i*TSIL_A_ (ma , Q2);
    
    Ac = -i*TSIL_A_ (mc , Q2);
    
    Aw = -i*TSIL_A_ (mw , Q2);
    
    Az = -i*TSIL_A_ (mz , Q2);
    
    Bac = i*TSIL_B_ (ma, mc, s, Q2);
    
    Bcw = i*TSIL_B_ (mc, mw, s, Q2);
    
    Bcz = i*TSIL_B_ (mc, mz, s, Q2);
    
    Bwc = Bcw;
    Bzc = Bcz;
    Bca = Bac;
    
    dBwc = i*TSIL_dBds_(mw,mc,s,Q2);
    dBzc = i*TSIL_dBds_(mz,mc,s,Q2);
    dBca = i*TSIL_dBds_(mc,ma,s,Q2);
    
  }
  
  int init(Data data)
  {
    /*mw= data.M_w, mz = data.M_z ,ma = data.M_a, mc = data.M_chi, g2=data.g1;
     tW = acos(mw/mz);
     sw2 = TSIL_POW(sin(tW),2),S2TW=TSIL_POW(sin(tW),2), cw2 = TSIL_POW(cos(tW),2);
     sw = sin(tW), cw = cos(tW);
     C = TSIL_POW(g2,2)/(16.0L*PI*PI);
     i=Power(-1,0.5);
     Pi=PI;
     TSIL_REAL p = data.M_chi, Q2 =data.Q;*/
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
  
  
  void Supplements::add_derivatives(Data &data)
  {
    TSIL_REAL p = 1.0,M=1.0;//data.M_chi, M = data.M_chi;
    init(data);
    TSIL_COMPLEXCPP SE_0 = (SigmaM_0(M) + M*SigmaK_0(M,p) )*( SigmaK_0(M,p) + 2.0L*TSIL_POW(M,2)*d_SigmaK_0(M,p)+2.0L*M*d_SigmaM_0(M));
    TSIL_COMPLEXCPP SE_1 =(SigmaM_1(M) + M*SigmaK_1(M,p) )*( SigmaK_1(M,p) + 2.0L*TSIL_POW(M,2)*d_SigmaK_1(M,p)+2.0L*M*d_SigmaM_1(M));
    //data.SE_1 = data.SE_1+real(SE_0);
    //data.SE_2 = data.SE_2+real(SE_1);
    
    
    cout << "SE derivative 0 = " << SE_0<<endl;;
    cout << "SE derivative 1 = " << SE_1<<endl;;
    
  }
}

