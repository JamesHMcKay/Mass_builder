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
	// electron mass
	long double M_e = 5.10998928E-04;
	
	// number of components in multiplet
	int components = 3;
	
	// Fermi constant
	long double G_F=1.16637876E-05;

  long double f_pi=0.13041;
  long double V_ud=0.97425;
	
	
	
	public:
	Decays (){}  // defualt constructor
	Decays(Data data) : data(data) {};
	
	
	long double calc_lifetime(long double deltam);
	
	long double pion_channel(long double deltam);
	
	long double muon_channel(long double deltam);
	
	long double electron_channel(long double deltam);
	
	long double pion_decay();
};

#endif
