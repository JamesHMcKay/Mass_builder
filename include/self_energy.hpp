#ifndef SELF_ENERGY_H
#define SELF_ENERGY_H
#include "data.hpp"
#include <complex>
using namespace std;
typedef complex<double> dcomp;
class Self_energy
{
  private:
  public:
  Data data;
  Self_energy (){}
  void run_tsil(Data &data);
};
namespace tsil
{
  class Integrals
  {
    public:
    std::complex<long double> Aa, Ac, Ah, Aw, Az ;
    std::complex<long double> Baw, Bca, Bcc, Bcw, Bcz ;
    std::complex<long double> Bwh, Bww, Bwz, Bzh ;
    std::complex<long double> Bwa;
    std::complex<long double> Bac;
    std::complex<long double> Bwc;
    std::complex<long double> Bzc;
    std::complex<long double> Bhw;
    std::complex<long double> Bzw;
    std::complex<long double> Bhz;
    std::complex<long double>  i;
    double  MChi,  MChi2 ,  ma,  ma2 ,  mw,  mw2 ,  mz,  mz2 ,  mh, mh2 ;

    double p,Pi;
  double cw, cw2, sw, g2, sw2, g1, STW, CTW, S2TW, C2TW, Lam, v, d1, dg2 ;

    Integrals (){}
    void DoTSIL(Data data);
  };
}
namespace F1_0
{
 
void  SE_F1(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace F2_0
{
 
void  SE_F2(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace V1_0
{
 
void  SE_V1(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace V2_0
{
 
void  SE_V2(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace V2_1
{
 
void  SE_V2(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace V3_0
{
 
void  SE_V3(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
namespace V3_1
{
 
void  SE_V3(Data data, tsil::Integrals integral);

double SE_1();double SE_2();
}
#endif
