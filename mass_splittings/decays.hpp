#ifndef DECAYS_H
#define DECAYS_H

#include "data.hpp"

using namespace std;

 
#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

class Decays
{
	private:
	Data data;
	int method;
	
	long double Pi = PI;
	
	long double deltam;
	
	// muon mass
	long double M_u = 1.05658372E-01;
	// pion mass
	long double M_pi = 0.139570L;
	// charged kaon mass
	long double M_ka = 0.493677L;
	// electron mass
	long double M_e = 5.10998928E-04;
	
	// see https://link.springer.com/content/pdf/10.1007%2FJHEP07%282015%29074.pdf
	int j; // isospin I = (-j,-j+1, . . . , j-1, j)
	int I; // T3 eigenvalue
	int T;
	
	// Fermi constant
	long double G_F = 1.16637876E-05;

  long double f_pi = 0.1302; // in GeV
  long double f_k = 0.1556; // in GeV
  long double V_ud = 0.97417;
  
  long double V_us =  0.2248;
	
	
	
	public:
	Decays (){}  // defualt constructor
	Decays(Data data, int j) : data(data), j(j) {};
	
	
	long double calc_lifetime(long double deltam, int I = 0);
	
	void calc_lifetime(long double deltam, ofstream &out_file, int I = 0);
	
	long double pion_channel(long double deltam);
	
	long double kaon_channel(long double deltam);
	
	long double muon_channel(long double deltam);
	
	long double electron_channel(long double deltam);
	
};

#endif
