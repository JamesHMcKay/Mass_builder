// calculate the lifetime of the charged component as a function of the mass splitting
// these functions call the decay widths, add them together and invert the result to give a lifetime

#include "decays.hpp"

using namespace std;

long double Decays::calc_lifetime(long double deltam)
{
	long double gamma = pion_channel(deltam) + electron_channel(deltam) + muon_channel(deltam);
	long double lifetime = 1.0 / gamma;

	return lifetime;
}

void Decays::calc_lifetime(long double deltam, ofstream &out_file)
{
	long double gamma = pion_channel(deltam) + electron_channel(deltam) + muon_channel(deltam) + kaon_channel(deltam);
	out_file << " " << deltam;
	out_file << " " << 1.0 / gamma;
	
	out_file << " " << pion_channel(deltam)/gamma;
	out_file << " " << muon_channel(deltam)/gamma;
	out_file << " " << electron_channel(deltam)/gamma;
	out_file << " " << kaon_channel(deltam)/gamma << endl;
}



// decay channels go here
// most have two options for how to calculate for the sake of comparison for different methods
// use the method flag in the data strcuture to control this choice
// each decay width function checks if mass difference is sufficient, rather do this in the
// function itself rather than in the call above so we can then use these widths independently

// subroutine used in both pion_channel function
long double Decays::pion_decay()  // taken from http://www2.ph.ed.ac.uk/~playfer/PPlect12.pdf
{
	long double A = pow(G_F * V_ud * f_pi* M_u,2) * M_pi/(8.0*Pi);

	return A*pow(1.0L-pow(M_u/M_pi,2),2)/6.582119e-25;
}

// muon decay, from arXiv:0903.3381
long double Decays::muon_channel(long double deltam)
{
	
	if (deltam > M_u)
	{
		return 0.12 * electron_channel(deltam);   
	}
	else
	{
		return 0;
	}

}


// kaon decay
long double Decays::kaon_channel(long double deltam)
{
	
	if (deltam > M_ka)
	{
	  long double A = (pow(components,2)-1.0L)*pow(G_F * f_pi * V_us ,2)*pow(deltam,3)/(4.0*Pi);
	  return A * pow( 1.0L - pow(M_ka/deltam,2),0.5)/6.582119e-25;
	}
	else
	{
	  return 0;
	}
	
}


// pion decay
long double Decays::pion_channel(long double deltam)
{
	
	if (deltam > M_pi)
	{
	  long double A = (pow(components,2)-1.0L)*pow(G_F * f_pi * V_ud ,2)*pow(deltam,3)/(4.0*Pi);
	  return A * pow( 1.0L - pow(M_pi/deltam,2),0.5)/6.582119e-25;
	}
	else
	{
	  return 0;
	}
	
}

// from arXiv:0903.3381
long double Decays::electron_channel(long double deltam)
{
	
	if (deltam > M_e)
	{
	  long double x = pow(M_e/deltam,2);
	  long double corrections;
	  
	  corrections=(1.0L-8.0L*x-12*pow(x,2)*log(x)+8*pow(x,3)-pow(x,4));
	  
	  long double A=(pow(components,2)-1.0L)*pow(G_F,2)/(60.0L*pow(Pi,3));
	  
	  return corrections*A*pow(deltam,5)/6.582119e-25;
	}
	else
	{
		return 0;
	}

}
