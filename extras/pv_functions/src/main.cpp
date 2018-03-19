// Simple program for calculating the MDM mass splitting using analytic expressions for the Passarino-Veltman functions from hep-ph/9606211

#include "pv.hpp"
#include "mdm.hpp"
#include "figures.hpp"
#include "data.hpp"
#include "decays.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>

using namespace std;

/*
// do some consistency checks for the PV functions
void consistency_check()
{
long double mz=91.1876, M=1e8, mw= 80.385;
long double Q=100.0,g=0.65,cw2=0.77,sw2=0.23;

cout<< "performing consistency checks " << endl;
B_0 B0(Q);
B_1 B1(Q);

// calculate mass splitting from difference of self energies (simplified form)
dcomp delta_m =-(2.0L*B0(M,M,mw)+B1(M,M,mw))+ cw2 * ( 2.0L*B0(M,M,mz)+B1(M,M,mz)) + sw2 * (2.0L*B0(M,M,0)+B1(M,M,0));
delta_m=-(M*pow(g,2))*(1/(8*pow(3.14159,2)))*delta_m;
// calculate the high M>>m_z,w limit
long double theory=(1/(8*3.14159))*(pow(g,2))*(mw-cw2*mz);

cout << "delta m = " << real(delta_m) << " compared with approximation = " << theory << endl;

// calculate a simple difference and compare with high M limit
cout << "B0(M,M,mz)-B0(M,M,mw) = " << real(M*(B0(M,M,mz)-B0(M,M,mw))) << " compared with approximation = " << 3.14159*(mz-mw) << endl;
cout<< "Large M, B0(M,M,M) = " << real(B0(M,M,M)) << " compared with approximation =  " << 3.14159/pow(3,0.5)-2 << endl;
}
*/

int main(int argc, char* argv[])
{

Data data(argc,argv);

//consistency_check();  // uncomment to compare the PV functions to limiting cases, only the last check
// isn't quite right because there is ambiguity in the definition of B when we aren't considering differences, not worried about this when considering
// mass splittings though

cout << " ------------------------------------------------------------------------------------------- " << endl;
cout << " ---------- calculating minimal dark matter pole masses and decay lifetimes ---------------- " << endl;
cout << " ------------------------------------------------------------------------------------------- " << endl;

MDM mdm(data);
//Two_loop two_loop(data);
//cout << "partial two-loop from Yamada 2009 = " << two_loop.partial_two_loop()<< endl;


//two_loop.init_tsil(data.M_chi);


//cout << "subset of two loop diagrams = " << two_loop.get_one_loop_plus_0(10000.0,10000.0)-two_loop.get_one_loop_plus_1(10000.0,10000.0) << endl;


//mdm.set_pole_mass(); // no argument uses default setting to calculate tree-level mass using iterative pole mass

data=mdm.data; // reset the data class to use the new value for M_chi tree-level mass

// print out some stuff using the MDM class

/*
cout << "the following are computed for M_tree = " << data.M_chi << endl;
cout << "pole mass (Mn,Mc) using iteration = " << mdm.calculate_pole_mass_n() << " , " << mdm.calculate_pole_mass_c() << endl;
cout << "pole mass (Mn,Mc) explicit method = " << mdm.calculate_pole_mass_n_simple() << " , " << mdm.calculate_pole_mass_c_simple() << endl;
cout << "pole mass (Mn,Mc) FeynCalc method = " << mdm.calculate_pole_mass_n_alt() << " , " << mdm.calculate_pole_mass_c_alt() << endl;
cout<< "mass splitting (using iterative pole mass) " << mdm.calculate_mass_splitting()*1000 << " MeV" << endl;
double delta_m=real(mdm.Sigma_n(data.M_chi)-mdm.Sigma_c(data.M_chi));
cout << "mass splitting from simple pole masses = " << delta_m*1000 << " MeV" << endl;
cout<< "mass splitting between doubly charged and neutral components " << (mdm.calculate_pole_mass_cc_simple()-mdm.calculate_pole_mass_n_simple())*1000 << " MeV"<< endl;


cout<< "alternative mass splitting " << (mdm.calculate_pole_mass_c_alt()-mdm.calculate_pole_mass_n_alt())*1000 << " MeV"<< endl;


cout<< "alternative mass splitting simple " << real(mdm.SE_F5(data.M_chi)-mdm.SE_F6(data.M_chi))*1000 << " MeV"<< endl;
*/

/*
Decays decays(data,mdm);

cout << "pion decay lifetime " << (1.0/decays.pion_decay())*1e9 << " ns" << endl;

cout << "simple decay lifetime is = " << decays.lifetime_simple()*3e11 <<" mm/c " << decays.lifetime_simple()*1e9 <<" ns" <<endl;
cout << "decay lifetime is = " << decays.lifetime()*3e11<<" mm/c " <<  decays.lifetime()*1e9<<" ns"  << endl;


mdm.set_pole_mass(1); // calculate tree-level mass using explicit pole mass expression

data=mdm.data; // reset the data class to use the new value for M_chi tree-level mass

cout << " ----------------------------------------------- " << endl;

// print out some stuff using the MDM class

cout << "the following are computed for M_tree = " << data.M_chi << endl;
cout << "pole mass (Mn,Mc) using iteration = " << mdm.calculate_pole_mass_n() << " , " << mdm.calculate_pole_mass_c() << endl;
cout << "pole mass (Mn,Mc) explicit method = " << mdm.calculate_pole_mass_n_simple() << " , " << mdm.calculate_pole_mass_c_simple() << endl;
cout<< "mass splitting (using iterative pole mass) " << mdm.calculate_mass_splitting()*1000 << " MeV" << endl;
delta_m=real(mdm.Sigma_n(data.M_chi)-mdm.Sigma_c(data.M_chi));
cout << "mass splitting from simple pole masses = " << delta_m*1000 << " MeV" << endl;
cout<< "mass splitting between doubly charged and neutral components " << (mdm.calculate_pole_mass_cc_simple()-mdm.calculate_pole_mass_n_simple())*1000 << " MeV"<< endl;
Decays decays2(data,mdm);

cout << "pion decay lifetime " << (1.0/decays2.pion_decay())*1e9 << " ns" << endl;

cout << "simple decay lifetime is = " << decays2.lifetime_simple()*3e11 <<" mm/c " << decays2.lifetime_simple()*1e9 <<" ns" <<endl;
cout << "decay lifetime is = " << decays2.lifetime()*3e11<<" mm/c " <<  decays2.lifetime()*1e9<<" ns"  << endl;

//PV pv(100.0);
//pv.plot_B0(1,100.0,2000.0,200.0,300.0); // ( int : argument to vary, double : lower limit, double : upper limit, (fix remaining two arguments))
//pv.plot_B1(1,100.0,2000.0,1000,90);

Data data2(argc,argv);
Figures fig(data2);*/
//fig.plot_mass_splitting_unc();
//fig.plot_mass_unc(); // be careful if you run a function before this which changes the data structure, as then it will calculate for M not equal to your original input

//fig2.plot_two_loop();

//fig.plot_approx();
//fig.plot_derivatives();
//fig.plot_limits();
//fig.plot_pole_masses_n();
//fig.plot_pole_masses_c();
//fig.plot_SigmaK();

Figures fig(data);
//fig.plot_decay();
//fig.plot_mass_splittings();
fig.plot_limits();


}
