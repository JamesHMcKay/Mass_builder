#ifndef INTEGRATE_H
#define INTEGRATE_H


#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;


struct B0_integrand {
B0_integrand(long double p,long double m_1,long double m_2,long double Q, long double epsilon) : p(p),m_1(m_1),m_2(m_2), Q(Q), epsilon(epsilon) {}
long double operator()(double x) const {
dcomp i;
i=-1;
i=sqrt(i);
return real(log( ( (1.0L-x)*pow(m_1,2)+x*pow(m_2,2)-x*(1.0L-x)*pow(p,2)-i*epsilon)/ pow(Q,2)));
}

private:
long double p,m_1,m_2;
long double Q,epsilon;
};



struct B0_der_integrand {
B0_der_integrand(long double p,long double m_1,long double m_2,long double Q, long double epsilon) : p(p),m_1(m_1),m_2(m_2), Q(Q), epsilon(epsilon) {}
long double operator()(long double x) const {
dcomp i;
i=-1;
i=sqrt(i);
return -x*(1-x)/ real( ( (1.0L-x)*pow(m_1,2)+x*pow(m_2,2)-x*(1.0L-x)*pow(p,2)-i*epsilon));
}

private:
long double p,m_1,m_2;
long double Q,epsilon;
};


struct Quadrature{
//abstract base class
int n; // current level or refinement
virtual double next() = 0 ;
};

template<class T>
struct Trapzd : Quadrature {
T &func;
long double a, b , s;

Trapzd() {}; // constructor
Trapzd(T &funcc, const double aa, const double bb) :
func ( funcc) , a(aa) , b(bb) {n=0;}

double next()
{
  double x, tnm, sum , del;
  int it,j;
  n++;
  if (n==1)
  {
    return (s=0.5*(b-a)*(func(a)+func(b)));
  }
  else
  {
    for (it=1, j=1; j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0, j=0;j<it;j++,x+=del) sum +=func(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
  
}

};


template<class T>
double qtrap(T &func, const double a, const double b, const double eps=1.0e-10)
{
const int JMAX=30;
double s,olds=0.0;
Trapzd<T> t(func,a,b);
for (int j=0;j<JMAX;j++) {
    s=t.next();
    if (j > 5)
if (abs(s-olds) < eps*abs(olds) ||
(s == 0.0 && olds == 0.0)) return s;
olds=s; }
cout<< "too many steps " << endl;//throw("Too many steps in routine qtrap");
}






#endif