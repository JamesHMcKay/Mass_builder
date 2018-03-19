#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>


#include <fstream>

#include "pv.hpp"
#include "integrate.hpp"


// contains definitions of the two point PV functions B_0(p,m_1,m_2) and B_1(p,m_1,m_2) obtained in analytic form from hep-ph/9606211
// and routines for plotting these functions (with respect to either p, m_1, m_2, with the other variables held constant), see below
// we neglect the divergent term as this is subtracted in minimal subtraction, we don't need it here for just calculating
// the pole masses


using namespace std;

dcomp PV::fB(dcomp x)
{
return log(1.0L-x)-x*log(1.0L-(1.0L/x))-1.0L;
}

dcomp PV::B_0()
{
  dcomp i;
  i=-1;
  i=sqrt(i);
  dcomp s=pow(p,2)-pow(m_2,2)+pow(m_1,2);

  dcomp xp=(s+pow(  pow(s,2) - 4.0L * pow (p,2) * (pow(m_1,2)-i * epsilon) ,0.5)) / (2.0L * pow(p,2) );

  dcomp xm=(s-pow(  pow(s,2) - 4.0L * pow (p,2) * (pow(m_1,2)-i * epsilon) ,0.5)) / (2.0L * pow(p,2) );
  //return -(fB(xp) + fB(xm))- 2.0L*log(p/Q);
  return -(fB(xp) + fB(xm))- 2.0L*log(p) + 2.0L*log(Q);
}

dcomp PV::A_0 (long double m)
{

long double result;
if (m==0){result=0;}
else
{
  result = pow(m,2)*(1.0L - 2.0L*log(m)+ 2.0L*log(Q));
}


return result;
}

/*
dcomp PV::A_0 ()
{

long double result;
if (m_1==0){result=0;}
else
{
result = pow(m_1,2)*(1.0L - 2.0*log(m_1/Q));
}
return result;
}
*/
dcomp PV::B_1 ()
{

  //cout << " A_0(m_2) " << A_0(m_2) << " A_0(m_1) " << A_0(m_1) << " pow(p,2)*B_0() " << pow(p,2)*B_0() << " pow(m2,2)*B_0() " << pow(m_2,2)*B_0() << " pow(m1,2)*B_0() " << pow(m_1,2)*B_0() << endl;
  //cout << "result = " <<(A_0(m_2)-A_0(m_1)+(pow(p,2)+pow(m_1,2)-pow(m_2,2))*B_0()) << endl;
  return -sign*(0.5L/pow(p,2))*(A_0(m_2)-A_0(m_1)+(pow(p,2)+pow(m_1,2)-pow(m_2,2))*B_0());

}


dcomp PV::B_22()
{
dcomp result;


result = 0.5L*(A_0(m_1)+A_0(m_2)) + ( pow(m_1,2)+pow(m_2,2) - 0.5L * pow(p,2) )*B_0()+( (pow(m_2,2)-pow(m_1,2))/(2.0L*pow(p,2)) ) * ( A_0(m_2)-A_0(m_1)-(pow(m_2,2)-pow(m_1,2))*B_0() )+pow(m_1,2)+pow(m_2,2)-0.33333333L*pow(p,2);

result = result / 6.0L;
return result; // slow with this precision, but seems to be required for large values of the B0 arguments
}

dcomp PV::B_22_tilda()
{
dcomp result;


result = B_22()-0.25L*A_0(m_1)-0.25L*A_0(m_2);

return result; // slow with this precision, but seems to be required for large values of the B0 arguments
}






///////////////////////////////////////////
////// hard coded limiting functions, not yet finished //////
///////////////////////////////////////////

//long double PV::B_0(long double M,long double delta,long double m)
//{
//// for B0 of the form B(M,delta*M,m)
//
//long double result;
//
//long double m2=pow(m,2), M2=pow(M,2), d2=pow(delta,2);
//
//long double A = log(1-d2/M- ( -1+d2) );
//long double B = log(1-d2/M + pow( 1- 2*d2+pow(d2,4), 0.5));
//
//long double S= pow( 1- 2*d2+ pow(d2,2),0.5);
//
//
//result += 0;//-2+log(4);
//result += ((m2-M2-M*d2+S)/(2*M*d2)) * A;
//result += ((-m2+M2+M*d2+S)/(2*M*d2)) * B;
//
//result +=
//
//return
//
////}










////////////////////////////////////////////
///////////////  Figures ///////////////////
////////////////////////////////////////////

// first argument -- variable to vary over  1:p, 2:m_1, 3:m_2
// second and third, upper and lower bound
// a and b are the fixed values for the fixed arguments, in numerical order

void PV::plot_B0(int n, long double lower, long double upper,long double a, long double b)
{
p=a,m_1=a,m_2=p;
int p_step=0,m_1_step=0,m_2_step=0;

if (n==1)
{
p=lower,m_1=a,m_2=b;
p_step=1,m_1_step=0,m_2_step=0;
}
else if (n==2)
{
m_1=lower,p=a,m_2=b;
p_step=0,m_1_step=1,m_2_step=0;
}
else
{
p=a,m_1=a,m_2=p;
p_step=0,m_1_step=0,m_2_step=1;
}

//dcomp B0_f=B_0();
//cout<< "steps (m1,m2,p) = " << m_1_step << " " << m_2_step << " " << p_step << endl;
// plot B0
ofstream myfile;
myfile.open ("../Figures/data/B0.txt");
long double step=((upper-lower)/float(1000));
for (int i=0;i<1000;i++)
{
p=p+step*p_step;
m_1=m_1+m_1_step*step;
m_2=m_2+m_2_step*step;
//cout<< "(m1,m2,p) = " << m_1 << " " << m_2 << " " << p << endl;

myfile << step*float(i)+lower << " " << real(B_0()) << endl;

}
myfile.close();
system("python ../Figures/B0.py");


}



void PV::plot_B0_r(long double r,long double lower, long double upper, long double m)
{

ofstream myfile;
myfile.open ("../Figures/data/B0_r.txt");
long double step=((upper-lower)/float(1000));
m_2=m;
m_1=lower;
long double m_1_step;
for (int i=0;i<1000;i++)
{
m_1=m_1+step;
p=m_1;
long double B0=real(B_0());
p=r*p;
long double B0r=real(B_0());

myfile << step*float(i)+lower << " " << B0 << " " << B0r << endl;
}
myfile.close();
system("python ../Figures/B0_r.py");

}



void PV::plot_B1(int n, long double lower, long double upper,long double a, long double b)
{
p=a,m_1=a,m_2=p;
int p_step=0,m_1_step=0,m_2_step=0;

if (n==1)
{
p=lower,m_1=a,m_2=b;
p_step=1,m_1_step=0,m_2_step=0;
}
else if (n==2)
{
m_1=lower,p=a,m_2=b;
p_step=0,m_1_step=1,m_2_step=0;
}
else
{
p=a,m_1=a,m_2=p;
p_step=0,m_1_step=0,m_2_step=1;
}

//dcomp B0_f=B_0();
//cout<< "steps (m1,m2,p) = " << m_1_step << " " << m_2_step << " " << p_step << endl;
// plot B0
ofstream myfile;
myfile.open ("../Figures/data/B1.txt");
long double step=((upper-lower)/float(1000));
for (int i=0;i<1000;i++)
{
p=p+step*p_step;
m_1=m_1+m_1_step*step;
m_2=m_2+m_2_step*step;
//cout<< "(m1,m2,p) = " << m_1 << " " << m_2 << " " << p << endl;

myfile << step*float(i)+lower << " " << real(B_1()) << endl;

}
myfile.close();
system("python ../Figures/B1.py");


}
