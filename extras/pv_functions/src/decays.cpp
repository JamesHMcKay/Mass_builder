#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include "mdm.hpp"
#include "pv.hpp"
#include "decays.hpp"

using namespace std;


// calculate the lifetime of the charged component as a function of the mass splitting
// these functions call the decay widths, add them together and invert the result to give a lifetime


long double Decays::calc_lifetime(long double deltam)
{
long double gamma = pion_channel(deltam)+electron_channel(deltam)+muon_channel(deltam);
long double lifetime = 1.0 / gamma;
return lifetime;
}

// use the pole mass from the implicit iterative expression

long double Decays::lifetime()
{
long double MFc = mdm.calculate_pole_mass_c();
long double MFn = mdm.calculate_pole_mass_n();

return calc_lifetime(MFc-MFn);
}
// use the pole mass from the simple explicit expression
long double Decays::lifetime_simple()
{
long double MFc = mdm.calculate_pole_mass_c_simple();
long double MFn = mdm.calculate_pole_mass_n_simple();
return calc_lifetime(MFc-MFn);
}


long double Decays::lifetime2()
{
long double MFcc = mdm.calculate_pole_mass_cc();
long double MFn = mdm.calculate_pole_mass_c();

return calc_lifetime(MFcc-MFn);
}
// use the pole mass from the simple explicit expression
long double Decays::lifetime_simple2()
{
long double MFcc = mdm.calculate_pole_mass_cc_simple();
long double MFn = mdm.calculate_pole_mass_c_simple();
return calc_lifetime(MFcc-MFn);
}




// decay channels go here
// most have two options for how to calculate for the sake of comparison for different methods
// use the method flag in the data strcuture to control this choice
// each decay width function checks if mass difference is sufficient, rather do this in the
// function itself rather than in the call above so we can then use these widths independently


// subroutine used in both pion_channel function
long double Decays::pion_decay()  // taken from http://www2.ph.ed.ac.uk/~playfer/PPlect12.pdf
{


long double A = pow(data.G_F*data.V_ud*data.f_pi*data.M_u,2)*data.M_pi/(8.0*data.Pi);

return A*pow(1.0L-pow(data.M_u/data.M_pi,2),2)/6.582119e-25;


}



// muon
long double Decays::muon_channel(long double deltam)
{

if (deltam>data.M_u)
{
return 0.12*electron_channel(deltam);   // from arXiv:0903.3381
}
else
{
return 0;
}

}



// pion
long double Decays::pion_channel(long double deltam)
{

if (deltam>data.M_pi)
{
  long double A=(pow(data.components,2)-1.0L)*pow(data.G_F*data.f_pi*data.V_ud,2)*pow(deltam,3)/(4.0*data.Pi);
  return A*pow( 1.0L - pow(data.M_pi/deltam,2),0.5)/6.582119e-25;
}
else
{
  return 0;
}
}



// electron
long double Decays::electron_channel(long double deltam)
{
if (deltam>data.M_e)
{
  long double x=pow(data.M_e/deltam,2);
  long double corrections;
  corrections=(1.0L-8.0L*x-12*pow(x,2)*log(x)+8*pow(x,3)-pow(x,4));
  
  long double A=(pow(data.components,2)-1.0L)*pow(data.G_F,2)/(60.0L*pow(data.Pi,3));
  return corrections*A*pow(deltam,5)/6.582119e-25;
}
else
{
return 0;
}

}




