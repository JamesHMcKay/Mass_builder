#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>


#include <fstream>

#include "mdm.hpp"
#include "pv.hpp"
#include "figures.hpp"

using namespace std;


// use this subroutine to make up some "theory" curve that we expect the
// mass splitting to follow, for now it is just a consistency check to
// show that the difference in self energies is due to the arguments used
// in the PV functions (change line long double M0 = M to either Mc or Mn (doesn't matter))




long double Figures::theory(long double M)
{
  //compute each pole mass
  long double Q=100;
  MDM mdm(data);
  
  long double mz=91.1876, mw= 80.385;
  long double g=0.65,cw2=0.77,sw2=0.23;
  
  B_0 B0(Q);
  B_1 B1(Q);
  long double Mn=mdm.calculate_pole_mass_n();
  long double Mc=mdm.calculate_pole_mass_c();//+0.1;//delta_m_approx;
  
  long double M0=M;
  
  
  dcomp sk1 =(4.0L*B1(Mn,M0,mw)+2.0L);
  sk1 =-(pow(g,2))*(1.0L/(16.0L*pow(3.14159,2)))*sk1;
  
  dcomp sm1 =(8.0L*B0(Mn,M0,mw)-4.0L);
  sm1 =-(M0*pow(g,2))*(1.0L/(16.0L*pow(3.14159,2)))*sm1;
  
  dcomp sk2 =(cw2*B1(Mc,M0,mz) + sw2 * B1(Mc,M0,0)+B1(Mc,M0,mw)+1.0L);
  sk2 = -(pow(g,2))*(1.0L/(8.0L*pow(3.14159,2)))*sk2;
  
  dcomp sm2 =(cw2*B0(Mc,M0,mz) + sw2 * B0(Mc,M0,0)+B0(Mc,M0,mw)-1.0L);
  sm2 = -(M0*pow(g,2))*(1.0L/(4.0L*pow(3.14159,2)))*sm2;
  
  dcomp diff=(Mn*sk1+sm1)-(Mc*sk2+sm2);
  
  //long double delta=M-Mn;
  
  //long double diff2=(pow(g,2)/(8*3.1415))*(mw-cw2*mz)*pow(M/(M+delta),2)*(1-delta/M);
  
  return real(diff);//*pow(M/(M+delta),2)*pow((1-delta/M),3);;//real(diff);
  
  
  
}


// plot MDM mass splittings
void Figures::plot_mass_splittings()
{
  dcomp delta_m_alt,delta_m_approx,delta_m_it;
  long double Q = 100.0L, n;
  long double M;
  
  B_0 B0(Q);
  B_1 B1(Q);
  
  ofstream myfile;
  myfile.open ("../Figures/data/mass_splittings.txt");
  
  long double max_M = 1e5;
  long double min_M = 100;
  
  long double logMax = log10(max_M);
  long double logMin = log10(min_M);
  int pts = 100;
  
  for (int i=0;i<pts+1;i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    M = pow(10.0L,n);
    
    data.M_chi=M;
    data.Q = M;
    MDM mdm(data);
    
    delta_m_approx = mdm.calculate_pole_mass_c_simple() - mdm.calculate_pole_mass_n_simple();
    
    delta_m_it = (mdm.calculate_mass_splitting());
    
    //long double DeltaM = theory(M);
    myfile << M << " " << real(delta_m_approx) << " " << real(delta_m_it) << " " << mdm.calculate_pole_mass_n() << " " << mdm.calculate_pole_mass_c() << endl;
    
  }
  myfile.close();
  
  system("python ../Figures/mass_splittings.py");
}

// plot difference pole mass and tree level mass as functions of tree level mass

void Figures::plot_pole_masses_n()
{
  long double M, n;
  
  ofstream myfile;
  myfile.open ("../Figures/data/pole_masses_n.txt");
  for (int i=0;i<100;i++)
  {
    n=(float(i)/float(100.0L))*5.0L;  //(float(i)/100)*2000+1;
    M=pow(10.0L,n)+0.1L;
    data.M_chi=M;
    MDM mdm(data);
    
    myfile << M << " " << mdm.calculate_pole_mass_n() << " " << mdm.calculate_pole_mass_n_simple() << endl;
  }
  myfile.close();
  system("python ../Figures/pole_masses_n.py");
}

void Figures::plot_pole_masses_c()
{
  long double M, n;
  
  ofstream myfile;
  myfile.open ("../Figures/data/pole_masses_c.txt");
  for (int i=0;i<100;i++)
  {
    n=(float(i)/float(100))*4;  //(float(i)/100)*2000+1;
    M=pow(10.0L,n)+0.1L;
    data.M_chi=M;
    MDM mdm(data);
    
    myfile << M << " " << mdm.calculate_pole_mass_c() << " " << mdm.calculate_pole_mass_c_simple() << endl;
  }
  myfile.close();
  system("python ../Figures/pole_masses_c.py");
}


void Figures::plot_SigmaK()
{
  long double M, n;
  
  ofstream myfile;
  myfile.open ("../Figures/data/Sigma_K0.txt");
  for (int i=0;i<100;i++)
  {
    n=(float(i)/float(100.0L))*6.0L;  //(float(i)/100)*2000+1;
    M=pow(10.0L,n)+0.1L;
    data.M_chi=M;
    MDM mdm(data);
    
    myfile << M << " " <<  real(mdm.Sigma_M0(M)) << endl;
  }
  myfile.close();
  system("python ../Figures/Sigma_K0.py");
}



void Figures::plot_limits()
{
  ofstream myfile;
  myfile.open ("../Figures/data/limits.txt");
  
  long double Q=1000,q,q1,q2,q3,q4,q5;
  B_0 B0(Q);
  B_1 B1(Q);
  long double ma = 100, mb = 200;
  long double n,M;
  for (int i=1;i<100;i++)
  {
    n=(float(i)/float(100))*8;  //(float(i)/100)*2000+1;
    M=pow(10,n);
    q1=real(M*(B0(M,M,ma)-B0(M,M,mb))/(3.14159L));
    q=0.9999;
    q2=real(M*(B0(M*q,M,ma)-B0(M*q,M,mb))/(3.14159L));
    q=0.99;
    q3=real(M*(B0(M*q,M,ma)-B0(M*q,M,mb))/(3.14159L));
    q=1.001;
    q4=real(M*(B0(M*q,M,ma)-B0(M*q,M,mb))/(3.14159L));
    q=1.0001;
    q5=real(M*(B0(M*q,M,ma)-B0(M*q,M,mb))/(3.14159L));
    
    
    myfile << M << " " << q1 << " " << q2 << " " << q3 << " " << q4 << " " << q5 << endl;
    
  }
  myfile.close();
  system("python ../Figures/limits.py");
}





void Figures::plot_derivatives()
{
  ofstream myfile;
  myfile.open ("../Figures/data/derivatives.txt");
  
  long double approx;
  long double M_pole;
  long double n,M;
  for (int i=1;i<100;i++)
  {
    n=(float(i)/float(100.0L))*3.0L;  //(float(i)/100)*2000+1;
    M=pow(10.0L,n);
    data.M_chi=M;
    MDM mdm(data);
    
    approx=mdm.Sigma_M0_der(M);
    
    
    myfile << M << " "  << approx  << " " << M_pole << endl;
    
  }
  myfile.close();
  system("python ../Figures/derivatives.py");
}





void Figures::plot_approx()
{
  ofstream myfile;
  myfile.open ("../Figures/data/approx.txt");
  
  long double approx,approx2,Sigma_M0,Sigma_K0,Sigma_Mpole, Sigma_Kpole;
  long double approxK,approxK2;
  long double M_pole;
  long double n,M;
  for (int i=1;i<100;i++)
  {
    n=(float(i)/float(100))*4;  //(float(i)/100)*2000+1;
    M=pow(10,n);
    data.M_chi=M;
    MDM mdm(data);
    M_pole=(mdm.calculate_pole_mass_n());
    
    Sigma_M0=real(mdm.Sigma_M0(M));
    Sigma_K0=real(mdm.Sigma_K0(M));
    Sigma_Mpole=real(mdm.Sigma_M0(M_pole));
    Sigma_Kpole=real(mdm.Sigma_K0(M_pole));
    
    approx=Sigma_M0+2.0L*M*(M_pole-M)*mdm.Sigma_M0_der(M);
    approxK=Sigma_K0+2.0L*M*(M_pole-M)*mdm.Sigma_K0_der(M);
    
    approx2=Sigma_M0-2.0L*M*(M*Sigma_K0+Sigma_M0)*mdm.Sigma_M0_der(M);
    approxK2=Sigma_K0-2.0L*M*(M*Sigma_K0+Sigma_M0)*mdm.Sigma_K0_der(M);
    
    myfile << M << "      " << Sigma_Mpole << "          " << approx << "         " << approx2 << "          " << Sigma_Kpole << "          " << approxK << "          " << approxK2 << endl;
    
  }
  myfile.close();
  system("python ../Figures/approx.py");
}


void Figures::plot_decay()
{
  long double M=1e4;
  long double n;
  
  ofstream myfile;
  myfile.open ("../Figures/data/decay.txt");
  
  ofstream myfile2;
  myfile2.open ("../Figures/data/decay2.txt");
  for (int i=0;i<1000;i++)
  {
    n=(float(i)/float(1000.0L))*5.0L;  //(float(i)/100)*2000+1;
    M=pow(10.0L,n);
    data.M_chi=M;
    MDM mdm(data);
    Decays decays(data,mdm);
    
    myfile << M << " " << decays.lifetime() << " " << decays.lifetime_simple() << endl;
    
    
    myfile2 << M << " " << decays.lifetime2() << " " << decays.lifetime_simple2() << endl;
  }
  myfile.close();
  myfile2.close();
  system("python ../Figures/decay.py");
}



void Figures::plot_mass_unc()
{
  ofstream myfile;
  myfile.open ("../Figures/data/mass_unc.txt");
  
  cout << "mass unc plot for input mass = " << data.M_chi << endl;
  for (int i=0;i<1000;i++)
  {
    long double Q=(float(i)/1000.0L)*400.0L+50.0L;
    data.Q=Q;
    MDM mdm(data);
    myfile << Q<< " " << mdm.calculate_pole_mass_n() << " " << mdm.calculate_pole_mass_n_simple() << endl;
    
  }
  myfile.close();
  system("python ../Figures/mass_unc.py");
  
}

void Figures::plot_mass_splitting_unc()
{
  ofstream myfile;
  myfile.open ("../Figures/data/mass_splitting_unc.txt");
  
  cout << "mass unc plot for input mass = " << data.M_chi << endl;
  for (int i=0;i<1000;i++)
  {
    long double Q=(float(i)/1000.0L)*400L+50.0L;
    data.Q=Q;
    
    MDM mdm(data);
    
    long double delta_m_approx=mdm.calculate_pole_mass_c_simple()-mdm.calculate_pole_mass_n_simple();
    
    long double delta_m_it=(mdm.calculate_mass_splitting());
    
    
    myfile << Q << " " << real(delta_m_approx) << " " << real(delta_m_it) << endl;
    
  }
  
  myfile.close();
  system("python ../Figures/mass_splitting_unc.py");
  
}
