#ifndef PV_H
#define PV_H

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>

#include <fstream>
using namespace std;
typedef complex<long double> dcomp;

class PV
{
  private:
  long double p,m_1,m_2,Q;
  long double epsilon=1e-5;
  long double sign=1; // sign of PV functions convention
  public:
  long double m_s;
  PV (){}  // defualt constructor
  PV (long double _Q){Q=_Q;}  // constructor if just making figures
  PV(long double _p, long double _m_1, long double _m_2,double _Q)
  {
      m_1=_m_1;
      m_2=_m_2;
      p=_p;
      Q=_Q;
   } //constructor
  
  
   dcomp B_0();
  
   dcomp B_22();
  
   dcomp B_22_tilda();
  
   dcomp A_0(long double m);
  
   dcomp A_0();

   dcomp fB(dcomp x);
  
   dcomp B_1();
  
   long double B_0_integral();
  
   long double B_0_der_integral();
  
   void plot_B0(int l2,long double l1,long double u1,long double u2,long double l3);
   void plot_B1(int l2,long double l1,long double u1,long double u2,long double l3);
  
   void plot_B0_r(long double r,long double l1,long double u1,long double a);
  
  
};


struct A_0 {
A_0(long double Q) : Q(Q) {}
dcomp operator()(long double p, long double m_1) const {
PV pv(p,m_1,0.0L,Q);
return pv.A_0();
}
private:
long double p,m_1,Q;
};


struct B_0 {
B_0(long double Q) : Q(Q) {}
dcomp operator()(long double p, long double m_1, long double m_2) const {
PV pv(p,m_1,m_2,Q);
return pv.B_0();
}
private:
long double p,m_1,m_2,Q;
};

struct B_1 {
B_1(long double Q) : Q(Q) {}
dcomp operator()(long double p, long double m_1, long double m_2) const {
PV pv(p,m_1,m_2,Q);
return pv.B_1();
}
private:
long double p,m_1,m_2,Q;
};

struct B_0_der_integral {
B_0_der_integral(long double Q) : Q(Q) {}
long double operator()(long double p, long double m_1,long  double m_2) const {
PV pv(p,m_1,m_2,Q);
return pv.B_0_der_integral();
}
private:
long double p,m_1,m_2,Q;
};


struct B_22 {
B_22(long double Q) : Q(Q) {}
dcomp operator()(long double p, long double m_1, long double m_2) const {
PV pv(p,m_1,m_2,Q);
return pv.B_22();
}
private:
long double p,m_1,m_2,Q;
};

struct B_22_tilda {
B_22_tilda(long double Q) : Q(Q) {}
dcomp operator()(long double p, long double m_1, long double m_2) const {
PV pv(p,m_1,m_2,Q);
return pv.B_22_tilda();
}
private:
long double p,m_1,m_2,Q;
};




#endif
