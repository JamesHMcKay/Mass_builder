/*
 Mass Builder
 
 -- EW_triplet.cpp --

 Contains functions for computing the EW triplet spectrum from FlexibleSUSY
 and additional interface to TSIL for derivatives of one-loop
 energies
  
 */

#include "ew_triplet_spectrum.hpp"

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

//inline double sqr(double a) { return a * a; }
#include "flexiblesusy/src/utils.h"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"

#include "flexiblesusy.hpp"
#define ALGORITHM1 Two_scale

#include "flexiblesusy/models/EW_triplet/EW_triplet_input_parameters.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_slha_io.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_spectrum_generator.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_two_scale_model.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_two_scale_model_slha.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_physical.hpp"
#include "flexiblesusy/models/EW_triplet/EW_triplet_info.hpp"

 
using namespace flexiblesusy;
using namespace softsusy;


using namespace std;

namespace extra_TSIL_interface_EW_triplet
{
#include TSIL_PATH
  using namespace std;
  TSIL_DATA bar;  
  TSIL_REAL Q2,Q;
  TSIL_REAL p;
  TSIL_COMPLEXCPP Log(TSIL_REAL a){complex<double> s(a/Q2,-0.000);return log(s);}
  TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}
  TSIL_COMPLEXCPP Power(TSIL_COMPLEXCPP a, int b){return pow(a,b);}
  TSIL_REAL Sin(TSIL_REAL a){return sin(a);}
  TSIL_REAL Cos(TSIL_REAL a){return cos(a);}
  TSIL_REAL Dot(TSIL_REAL a, TSIL_REAL b){return a*b;}
  TSIL_REAL Sqrt(TSIL_REAL a){return TSIL_POW(a,0.5);}
  TSIL_COMPLEXCPP Ae(TSIL_REAL a) { return TSIL_Aeps_(TSIL_POW(a,2),Q2);}
  TSIL_COMPLEXCPP Be(TSIL_REAL a, TSIL_REAL b) { return TSIL_Beps_(TSIL_POW(a,2),TSIL_POW(b,2), TSIL_POW(p,2), Q2);}
  int          init(Data data);
  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}
  TSIL_COMPLEXCPP operator+(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}
  TSIL_COMPLEXCPP operator+(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}
  TSIL_COMPLEXCPP operator-(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}
  TSIL_COMPLEXCPP operator-(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}
  TSIL_COMPLEXCPP operator/(TSIL_COMPLEXCPP a,double b){TSIL_COMPLEXCPP c=b;return a/c;}
  TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}
  
  std::complex<long double>  i;
  long double  MChi,  MChi2 ,  ma,  ma2 ,mt , mt2,  mw,  mw2 ,  mz,  mz2 ,  mh,  mh2, mf, mf2 ,  null, null2 ;
  
  double Pi;
  long double alpha,CTW, Ca, Cw1, Cw2, Cz1, Cz2, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, MassBuilderJEpsilon, STW, cw, cw2, d1Z, d1m, d2Z, d2m, dMWsq1, dMZsq1, dZAA1, dZAZ1, dZW1, dZZA1, dZZZ1, dg2, e, g, g1, g2, sw, sw2, v ;
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca,Bac,Bcz,Bcw;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,dBac,M;
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az;
  
  
  void DoTSIL(Data data)
  {
    alpha=data.alpha, cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   STW = data.STW,   CTW = data.CTW,   e = data.e,   d1Z = data.d1Z,   d1m = data.d1m,   d2Z = data.d2Z,   d2m = data.d2m,   dZW1 = data.dZW1,   dMWsq1 = data.dMWsq1,   dMZsq1 = data.dMZsq1,   dZAA1 = data.dZAA1,   dZZA1 = data.dZZA1,   dZAZ1 = data.dZAZ1,   dZZZ1 = data.dZZZ1,   dg2 = data.dg2,   Ca = data.Ca,   Cz1 = data.Cz1,   Cz2 = data.Cz2,   Cw1 = data.Cw1,   Cw2 = data.Cw2,   v = data.v,   MassBuilderJEpsilon = data.MassBuilderJEpsilon ;
    
    TSIL_REAL Q2 = pow(data.Q,2);
    TSIL_REAL s = pow(data.P,2);
    MChi = data.MChi, MChi2 = TSIL_POW(data.MChi, 2) ,   ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,   mw = data.mw, mw2 = TSIL_POW(data.mw, 2) ,   mz = data.mz, mz2 = TSIL_POW(data.mz, 2) ,   mh = data.mh, mh2 = TSIL_POW(data.mh, 2) ,   null = data.null, null2 = TSIL_POW(data.null, 2) ,   mf = data.mf, mf2 = TSIL_POW(data.mf, 2) ,   mt = data.mt , mt2 = TSIL_POW(data.mt, 2) ;
    
    dcomp ii=-1;ii=sqrt(ii);i=ii;
    Pi=PI;
    cw =  mw/mz;
    cw2 =  TSIL_POW(cw,2);
    sw =  TSIL_POW(1.-cw2,0.5);
    sw2 =  TSIL_POW(sw,2);
    STW =  sw;
    CTW =  cw;
    e =  TSIL_POW(4*Pi*alpha,0.5);
    
    TSIL_REAL a = -1., b = 1., c = 1.;
    
    Aa = a*TSIL_A_ (ma2 , Q2);
    
    Ac = a*TSIL_A_ (MChi2 , Q2);
    
    Aw = a*TSIL_A_ (mw2 , Q2);
    
    Az = a*TSIL_A_ (mz2 , Q2);
    
    Bac = b*TSIL_B_ (ma2, MChi2, s, Q2);
    
    Bcw = b*TSIL_B_ (MChi2, mw2, s, Q2);
    
    Bcz = b*TSIL_B_ (MChi2, mz2, s, Q2);
    
    Bwc = Bcw;
    Bzc = Bcz;
    Bca = Bac;
    
    dBwc = c*TSIL_dBds_(mw2,MChi2,s,Q2);
    dBzc = c*TSIL_dBds_(mz2,MChi2,s,Q2);
    dBca = c*TSIL_dBds_(ma2,MChi2,s,Q2);
    
    dBac = c*TSIL_dBds_(ma2,MChi2,s,Q2);
    
  }
  
  double add_derivatives(Data &data)
  {
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*MChi*(-(((Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
                                                  (Ac*Power(MChi,2) - Aw*Power(MChi,2) - Power(MChi,4) - 4*dBwc*Power(MChi,6) + Bcw*Power(MChi,2)*Power(mw,2) -
                                                   2*dBwc*Power(MChi,4)*Power(mw,2)))/Power(MChi,6)) +
                                               ((2 + 8*dBwc*Power(MChi,2) + 8*Power(cw,2)*dBzc*Power(MChi,2) -
                                                 (2*(Ac - Aw - Bcw*Power(MChi,2) + 2*dBwc*Power(MChi,4) + Bcw*Power(mw,2) - dBwc*Power(MChi,2)*Power(mw,2)))/Power(MChi,2) +
                                                 (Ac - Aw - Bcw*(2*Power(MChi,2) - Power(mw,2)))/Power(MChi,2) -
                                                 (2*Power(cw,2)*(Ac - Az - Bcz*Power(MChi,2) + 2*dBzc*Power(MChi,4) + Bcz*Power(mz,2) - dBzc*Power(MChi,2)*Power(mz,2)))/
                                                 Power(MChi,2) + (Power(cw,2)*(Ac - Az - Bcz*(2*Power(MChi,2) - Power(mz,2))))/Power(MChi,2) +
                                                 8*dBac*Power(MChi,2)*Power(sw,2) - (2*(-Aa + Ac + Bca*Power(ma,2) - Bca*Power(MChi,2) - dBac*Power(ma,2)*Power(MChi,2) +
                                                                                        2*dBac*Power(MChi,4))*Power(sw,2))/Power(MChi,2) +
                                                 ((-Aa + Ac - Bca*(-Power(ma,2) + 2*Power(MChi,2)))*Power(sw,2))/Power(MChi,2))*
                                                (Aw + Az*Power(cw,2) + 2*Power(MChi,2) - 2*Bcw*Power(MChi,2) - 2*Bcz*Power(cw,2)*Power(MChi,2) - Bcw*Power(mw,2) -
                                                 Bcz*Power(cw,2)*Power(mz,2) + Aa*Power(sw,2) - Bca*Power(ma,2)*Power(sw,2) - 2*Bca*Power(MChi,2)*Power(sw,2) -
                                                 Ac*(1 + Power(cw,2) + Power(sw,2))))/(4.*Power(MChi,2))))/(64.*Power(Pi,4)*Power(sw,4));
    
    return real(result);
  }
  
  
  double F11_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*(Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
                              (-Ac + Aw + Power(MChi,2) + 4*dBwc*Power(MChi,4) - Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2)))/
    (64.*Power(MChi,3)*Power(Pi,4)*Power(sw,4));
    
    return real(result);
    
  }
  
  
  
  double F12_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(e,4)*(Aw*Power(cw,2) + Power(MChi,2) + Power(cw,2)*Power(MChi,2) + 4*Power(cw,2)*dBwc*Power(MChi,4) + 4*dBzc*Power(MChi,4) -
                                          Bcw*Power(cw,2)*Power(mw,2) + 2*Power(cw,2)*dBwc*Power(MChi,2)*Power(mw,2) - Bcz*Power(mz,2) + 2*dBzc*Power(MChi,2)*Power(mz,2) +
                                          Aa*Power(cw,2)*Power(sw,2) - Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) - 2*Power(MChi,2)*Power(sw,2) +
                                          Power(cw,2)*Power(MChi,2)*Power(sw,2) + 2*Power(cw,2)*dBac*Power(ma,2)*Power(MChi,2)*Power(sw,2) +
                                          4*Power(cw,2)*dBac*Power(MChi,4)*Power(sw,2) - 8*dBzc*Power(MChi,4)*Power(sw,2) + 2*Bcz*Power(mz,2)*Power(sw,2) -
                                          4*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,2) + Power(MChi,2)*Power(sw,4) + 4*dBzc*Power(MChi,4)*Power(sw,4) -
                                          Bcz*Power(mz,2)*Power(sw,4) + 2*dBzc*Power(MChi,2)*Power(mz,2)*Power(sw,4) + Az*Power(-1 + Power(sw,2),2) -
                                          Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2))))*
                              (-(Aw*Power(cw,2)) - Power(MChi,2) + 2*Bcz*Power(MChi,2) - Power(cw,2)*Power(MChi,2) + 2*Bcw*Power(cw,2)*Power(MChi,2) +
                               Bcw*Power(cw,2)*Power(mw,2) + Bcz*Power(mz,2) - Aa*Power(cw,2)*Power(sw,2) + Bca*Power(cw,2)*Power(ma,2)*Power(sw,2) +
                               2*Power(MChi,2)*Power(sw,2) - 4*Bcz*Power(MChi,2)*Power(sw,2) - Power(cw,2)*Power(MChi,2)*Power(sw,2) +
                               2*Bca*Power(cw,2)*Power(MChi,2)*Power(sw,2) - 2*Bcz*Power(mz,2)*Power(sw,2) - Power(MChi,2)*Power(sw,4) +
                               2*Bcz*Power(MChi,2)*Power(sw,4) + Bcz*Power(mz,2)*Power(sw,4) - Az*Power(-1 + Power(sw,2),2) +
                               Ac*(Power(-1 + Power(sw,2),2) + Power(cw,2)*(1 + Power(sw,2)))))/(256.*Power(cw,4)*Power(MChi,3)*Power(Pi,4)*Power(sw,4));
    
    return real(result);
    
  }
  
  
  
  
  TSIL_COMPLEXCPP  gammagamma_Chi(Data data,double Q)
  {
    DoTSIL(data);
    
    p = Q;// TSIL_POW(data.Q,0.5);
    TSIL_REAL Q2 = pow(Q,2);
    
    TSIL_COMPLEXCPP AcMB = -i*TSIL_A_ (MChi2 ,  Q2);
    
    // evaluate as s = Q^2
    TSIL_COMPLEXCPP BccMB = i*TSIL_B_ (MChi2, MChi2, Q2, Q2);
    
    TSIL_COMPLEXCPP C0 = (4*Power(e,2)*(6*Power(MChi,2) - (-Power(p,2))))/9.;
    TSIL_COMPLEXCPP CAc =   Complex(0,2.6666666666666665)*Power(e,2) ;
    TSIL_COMPLEXCPP CBcc =   Complex(0,-1.3333333333333333)*Power(e,2)*(2*Power(MChi,2) + Power(p,2)) ;
    
    TSIL_COMPLEXCPP result = + C0  + AcMB * CAc + BccMB * CBcc;
    
    return -result/(16.0L*TSIL_POW(PI,2));
  }
  
}

double EW_triplet_spectrum::add_1loop_der(Data &data)
{
	
	return extra_TSIL_interface_EW_triplet::add_derivatives(data);
	
}


double EW_triplet_spectrum::iterative_ms_bar_mass(Data data, string particle)
{
  long double M_tree = data.MChi;
  long double M_pole = data.MChi;
  data.P = M_tree;
  long double diff = 1;
  long double precision = 1e-8;
  int iteration = 0;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    if (data.do_tsil_all)
    {
			M_pole = M_tree + data.SE_1[particle] + data.SE_2[particle];
    }
    else
    {
      M_pole = M_tree + data.SE_1[particle];
    }
    
    diff = abs(M_pole - data.P);
    data.P = M_pole;
    
    //cout<< "\r" << "M_pole - p = " << diff << " GeV";
    //std::cout << std::flush;
    
    iteration++;
  } while (diff > precision  && iteration < 500);
  
  //cout<< "\r" << "M_pole - p = " << diff << " GeV";
  //cout << "\n";
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
    return 0;
  }
    
  return M_pole;
}

// determine MSbar parameters
void EW_triplet_spectrum::compute_spectra_MB_1loop()
{
  Self_energy se;
  
  // RGE evolution
  double A = (19./(10.*Pi)); // EW_triplet
  //double A = (17./(30.*Pi)); // Wino model
  //double A = (7./(30.*Pi)); // SM
  double alpha_mz = data.alpha;
  double mu0 = data.mz;
  double mu = data.Q;//pow(data.Q,0.5);
  
  data.alpha = pow(   1.0/alpha_mz -  A * log( mu/mu0) , -1);
  
  data.do_tsil_all = false;
  
}

// determine MSbar parameters
void EW_triplet_spectrum::compute_spectra_MB_2loop()
{
  Self_energy se;
  
  // matching (?) of SM to Wino model
  data.alpha = data.alpha*( 1.0-real(extra_TSIL_interface_EW_triplet::gammagamma_Chi(data,data.Q)) /pow(data.Q,2) );
  
  // determine MS bar input parameters at Q
  // only recompute the 1-loop functions for efficiency
  data.do_tsil_all = false;
  
  data.P = data.mw;
  se.run_tsil(data);
  data.mw = pow( pow(data.mw,2) - real(data.SE_1["V3"]),0.5 );
  
  data.P = data.mz;
  se.run_tsil(data);
  data.mz = pow( pow(data.mz,2) - real(data.SE_1["V2"]) ,0.5);
  
  // need to iterate to determine MS bar mass for MChi to match equation (9) of Ibe et al.
  //data.MChi = iterative_ms_bar_mass(data, "F11_g1");
  
  data.P = data.MChi;
  data.do_tsil_all = true;
}


bool EW_triplet_spectrum::compute_spectra_flexiblesusy()
{
	Spectrum_generator_settings spectrum_generator_settings;
  
  QedQcd oneset;
  EW_triplet_input_parameters input;
  EW_triplet_spectrum_generator<Two_scale> spectrum_generator;
  
  EW_triplet_slha_io slha_io;
  
  const std::string slha_file="flexiblesusy/models/EW_triplet/LesHouches.in.EW_triplet";
  
  slha_io.read_from_file(slha_file);
  slha_io.fill(oneset);
  slha_io.fill(input);
  slha_io.fill(spectrum_generator_settings);
  
  spectrum_generator.set_settings(spectrum_generator_settings);
  spectrum_generator.set_parameter_output_scale(slha_io.get_parameter_output_scale());
  
  //oneset.setPoleMt(data.mt);
  
  oneset.toMz();
  
	input.QEWSB=data.Q;
	input.Qin=data.Q;
  input.HiggsIN = 0.5*pow(data.mh,2);
  input.YcIN = 0.5*data.MChi;
  
  
  spectrum_generator.run(oneset, input);
	
	std::ostringstream warnings;
	const Problems<EW_triplet_info::NUMBER_OF_PARTICLES>& problems
	= spectrum_generator.get_problems();
	const bool error = problems.have_problem();
	problems.print_warnings(warnings);
	if (error==1)
	{
		// check for errors
		std::ostringstream problems_str;
		problems.print_problems(problems_str);
		
		cout<< FORMAT_SPINFO(4,problems_str.str()) << endl;
		
		return error;
	}
	EW_triplet_slha<Two_scale> model(spectrum_generator.get_model());
	  
  // update data struct with computed spectrum
  
  // MS bar masses
  
  data.mw = model.get_MVWp();
  data.mz = model.get_MVZ();
  data.mh = model.get_Mhh();
  data.mt = model.get_MFu(2);
  data.v = model.get_v();
  
  double g1 = pow(3./5.,0.5)*model.get_g1();
  double g2 = model.get_g2();
  
  data.alpha = pow(g1*g2,2) / (4 * Pi * (g1*g1 + g2*g2));
  
  //data.sw = Sin(model.get_ZZ(0,0)); // not used
  
	data.SE_1["F11_g1"] = model.get_MFn_pole_slha() - data.MChi;
	data.SE_1["F12_g1"] = model.get_MFc_pole_slha() - data.MChi;
  
  return error;
}


void EW_triplet_spectrum::compute_tsil()
{
	Self_energy se;
	se.run_tsil(data);	
	
	// add derivatives of 1-loop self energies
	data.SE_2["F11_g1"] = data.SE_2["F11_g1"] +  extra_TSIL_interface_EW_triplet::F11_der(data);
	data.SE_2["F12_g1"] = data.SE_2["F12_g1"] +  extra_TSIL_interface_EW_triplet::F12_der(data);	

	if (!data.do_tsil_all)
	{
		data.SE_2["F11_g1"] = 0;
		data.SE_2["F12_g1"] = 0;
	}

}


void EW_triplet_spectrum::compute_tsil_iterative()
{	
	Self_energy se;
	se.run_tsil(data);
	
	double m_F11 = iterative_ms_bar_mass(data,"F11_g1");
	double m_F12 = iterative_ms_bar_mass(data,"F12_g1");
	
	if ( (m_F11 == 0 ) || (m_F12 == 0 ) )
	{
		data.SE_1["F11_g1"] = 0;
		data.SE_1["F12_g1"] = 0;
		data.SE_2["F11_g1"] = 0;
		data.SE_2["F12_g1"] = 0;
	}
	else if (data.do_tsil_all)
	{
		data.SE_2["F11_g1"] = m_F11 - data.MChi;
		data.SE_2["F12_g1"] = m_F12 - data.MChi;
	}
	else
	{
		data.SE_1["F11_g1"] = m_F11 - data.MChi;
		data.SE_1["F12_g1"] = m_F12 - data.MChi;
	}	

}



double EW_triplet_spectrum::get_deltam()
{
	return data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
}

double EW_triplet_spectrum::get_deltam_2loop()
{
	return data.SE_1["F12_g1"]+data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - data.SE_1["F11_g1"];
}

double EW_triplet_spectrum::get_charged_mass()
{
	return data.MChi + data.SE_1["F12_g1"] ;
}

double EW_triplet_spectrum::get_neutral_mass()
{
	return data.MChi + data.SE_1["F11_g1"] ;
}

