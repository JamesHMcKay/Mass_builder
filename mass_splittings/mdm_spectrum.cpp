/*
 Mass Builder
 
 -- mdm_spectrum.cpp --

 Contains functions for computing MDM spectrum from FlexibleSUSY
 and additional interface to TSIL for derivatives of one-loop
 energies
  
 */

#include "mdm_spectrum.hpp"

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

//inline double sqr(double a) { return a * a; }
#include "flexiblesusy/src/utils.h"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"

#include "flexiblesusy.hpp"
#define ALGORITHM1 Two_scale

#include "flexiblesusy/models/MDM/MDM_input_parameters.hpp"
#include "flexiblesusy/models/MDM/MDM_slha_io.hpp"
#include "flexiblesusy/models/MDM/MDM_spectrum_generator.hpp"
#include "flexiblesusy/models/MDM/MDM_two_scale_model.hpp"
#include "flexiblesusy/models/MDM/MDM_two_scale_model_slha.hpp"
#include "flexiblesusy/models/MDM/MDM_physical.hpp"
#include "flexiblesusy/models/MDM/MDM_info.hpp"

 
using namespace flexiblesusy;
using namespace softsusy;


using namespace std;

namespace extra_TSIL_interface_MDM
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
  long double  MChi,  MChi2 ,  ma,  ma2 ,mt , mt2,  mw,  mw2 ,  mz,  mz2 ;
  
  double Pi;
  long double alpha,CTW, STW, cw, cw2, e, g, g1, g2, sw, sw2, v ;
  
  TSIL_COMPLEXCPP Bwc, Bzc, Bca,Bac,Bcz,Bcw;
  
  TSIL_COMPLEXCPP dBwc,dBzc,dBca,dBac,M;
  
  TSIL_COMPLEXCPP Aa, Ac, Aw, Az;
  
  
  void DoTSIL(Data data)
  {
    v= data.v, alpha=data.alpha, cw = data.cw,   cw2 = data.cw2,   sw = data.sw,   sw2 = data.sw2,   STW = data.STW,   CTW = data.CTW,   e = data.e;
    
    TSIL_REAL Q2 = pow(data.Q,2);
    TSIL_REAL s = pow(data.P,2);
    MChi = data.MChi, MChi2 = TSIL_POW(data.MChi, 2) ,   ma = data.ma, ma2 = TSIL_POW(data.ma, 2) ,   mw = data.mw, mw2 = TSIL_POW(data.mw, 2) ,   mz = data.mz, mz2 = TSIL_POW(data.mz, 2);
    
    dcomp ii=-1;ii=sqrt(ii);i=ii;
    Pi=PI;
    
    cw =  mw/mz;
    cw2 =  TSIL_POW(cw,2);
    sw =  TSIL_POW(1.-cw2,0.5);
    sw2 =  TSIL_POW(sw,2);
    STW =  sw;
    CTW =  cw;
    e =  TSIL_POW(4*Pi*alpha,0.5);
    g1 =  e/cw;
    g2 =  e/sw;
    v =  2*mw/g2;
    
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
 
  
  double F5_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(g2,4)*(5*Aw + Az*Power(CTW,2) + 5*Power(MChi,2) + Power(CTW,2)*Power(MChi,2) + 20*dBwc*Power(MChi,4) + 
       4*Power(CTW,2)*dBzc*Power(MChi,4) - 5*Bcw*Power(mw,2) + 10*dBwc*Power(MChi,2)*Power(mw,2) - Bcz*Power(CTW,2)*Power(mz,2) + 
       2*Power(CTW,2)*dBzc*Power(MChi,2)*Power(mz,2) + Aa*Power(STW,2) - Bca*Power(ma,2)*Power(STW,2) + Power(MChi,2)*Power(STW,2) + 
       2*dBac*Power(ma,2)*Power(MChi,2)*Power(STW,2) + 4*dBac*Power(MChi,4)*Power(STW,2) - Ac*(5 + Power(CTW,2) + Power(STW,2)))*
     (-5*Aw - Az*Power(CTW,2) - 5*Power(MChi,2) + 10*Bcw*Power(MChi,2) - Power(CTW,2)*Power(MChi,2) + 2*Bcz*Power(CTW,2)*Power(MChi,2) + 
       5*Bcw*Power(mw,2) + Bcz*Power(CTW,2)*Power(mz,2) - Aa*Power(STW,2) + Bca*Power(ma,2)*Power(STW,2) - Power(MChi,2)*Power(STW,2) + 
       2*Bca*Power(MChi,2)*Power(STW,2) + Ac*(5 + Power(CTW,2) + Power(STW,2))))/(256.*Power(MChi,3)*Power(Pi,4));
    
    return real(result);
    
  }
  
  
  
  double F6_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (Power(g2,4)*(Ac - Aw + 2*Ac*Power(CTW,2) - 2*Az*Power(CTW,2) - Power(MChi,2) + 2*Bcw*Power(MChi,2) - 2*Power(CTW,2)*Power(MChi,2) + 
       4*Bcz*Power(CTW,2)*Power(MChi,2) + Bcw*Power(mw,2) + 2*Bcz*Power(CTW,2)*Power(mz,2) - 2*Aa*Power(STW,2) + 2*Ac*Power(STW,2) + 
       2*Bca*Power(ma,2)*Power(STW,2) - 2*Power(MChi,2)*Power(STW,2) + 4*Bca*Power(MChi,2)*Power(STW,2))*
     (Aw + 2*Az*Power(CTW,2) + Power(MChi,2) + 2*Power(CTW,2)*Power(MChi,2) + 4*dBwc*Power(MChi,4) + 8*Power(CTW,2)*dBzc*Power(MChi,4) - 
       Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2) - 2*Bcz*Power(CTW,2)*Power(mz,2) + 
       4*Power(CTW,2)*dBzc*Power(MChi,2)*Power(mz,2) + 2*Aa*Power(STW,2) - 2*Bca*Power(ma,2)*Power(STW,2) + 
       2*Power(MChi,2)*Power(STW,2) + 4*dBac*Power(ma,2)*Power(MChi,2)*Power(STW,2) + 8*dBac*Power(MChi,4)*Power(STW,2) - 
       Ac*(1 + 2*Power(CTW,2) + 2*Power(STW,2))))/(64.*Power(MChi,3)*Power(Pi,4));
    return real(result);
    
  }
  
  double F7_der(Data data)
  {
    
    DoTSIL(data);
    
    p = data.P;
    
    
    TSIL_COMPLEXCPP result = (9*Power(g2,4)*(Ac - Aw - Power(MChi,2) + 2*Bcw*Power(MChi,2) + Bcw*Power(mw,2))*
     (-Ac + Aw + Power(MChi,2) + 4*dBwc*Power(MChi,4) - Bcw*Power(mw,2) + 2*dBwc*Power(MChi,2)*Power(mw,2)))/
   (64.*Power(MChi,3)*Power(Pi,4));
    
    return real(result);
    
  }
   
  
  TSIL_COMPLEXCPP  gammagamma_Chi(Data data,double Q)
  {
    DoTSIL(data);
    
    p = Q;// TSIL_POW(data.Q,0.5);
    TSIL_REAL Q2 = pow(Q,2);
    
    // evaluate as s = Q^2
    TSIL_COMPLEXCPP AcMB = -i*TSIL_A_ (MChi2 ,  Q2);
    TSIL_COMPLEXCPP BccMB = i*TSIL_B_ (MChi2, MChi2, Q2, Q2);
	     
	  TSIL_COMPLEXCPP C20 = (4*Power(g2,2)*(6*Power(MChi,2) - Power(p,2))*Power(STW,2))/9.;
	  TSIL_COMPLEXCPP C2AC =   Complex(0,2.6666666666666665)*Power(g2,2)*Power(STW,2) ;
	  TSIL_COMPLEXCPP C2BCC =   Complex(0,-1.3333333333333333)*Power(g2,2)*Power(STW,2)*(2*Power(MChi,2) + Power(p,2)) ;
	  TSIL_COMPLEXCPP result =  + C20  + AcMB * C2AC + BccMB * C2BCC;
	  
	  // use this if using full quintuplet model
	  
	  TSIL_COMPLEXCPP C20_2 = (16*Power(g2,2)*(6*Power(MChi,2) - Power(p,2))*Power(STW,2))/9.;
		TSIL_COMPLEXCPP C2AC_2 =   Complex(0,10.666666666666666)*Power(g2,2)*Power(STW,2) ;
		TSIL_COMPLEXCPP C2BCC_2 =   Complex(0,-5.333333333333333)*Power(g2,2)*Power(STW,2)*(2*Power(MChi,2) + Power(p,2)) ;
		TSIL_COMPLEXCPP result2 =  + C20_2  + AcMB * C2AC_2 + BccMB * C2BCC_2;

    return -(result/*+result2*/)/(16.0L*TSIL_POW(PI,2));
  }
  
}

double MDM_spectrum::iterative_ms_bar_mass(Data data, string particle)
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
void MDM_spectrum::compute_spectra_MB_1loop()
{
  Self_energy se;
  
  // RGE evolution
  double A = (19./(10.*Pi)); // MSSM
  //double A = (17./(30.*Pi)); // Wino model
  //double A = (7./(30.*Pi)); // SM
  double alpha_mz = data.alpha;
  double mu0 = data.mz;
  double mu = data.Q;//pow(data.Q,0.5);
  
  data.alpha = pow(   1.0/alpha_mz -  A * log( mu/mu0) , -1);
  
  data.do_tsil_all = false;
  
}

/*
// determine MSbar parameters
void MDM_spectrum::compute_spectra_MB_2loop()
{
	
	cout << "------------- Computing MB spectrum -------------" << endl;
  Self_energy se;  
  // matching (?) of SM to Wino model
  data.alpha = data.alpha*( 1.0-real(extra_TSIL_interface_MDM::gammagamma_Chi(data,data.Q)) /pow(data.Q,2) );
  
  // determine MS bar input parameters at Q
  // only recompute the 1-loop functions for efficiency
  data.do_tsil_all = false;
  cout << "MB: mw , mz = " << data.mw << ", " << data.mz << endl;
  data.P = data.mw;
  se.run_tsil(data);
  data.mw = pow( pow(data.mw,2) - real(data.SE_1["V3"]),0.5 );

	double mz_pole = data.mz;
  
  data.P = 91.1876;
  se.run_tsil(data);
  data.mz = pow( pow(data.mz,2) - real(data.SE_1["V2"]) ,0.5);
  
  std::cout << "MB: MZ pole = " << mz_pole << " ms_bar = " << data.mz << std::endl;
  
  cout << "MB: self energy = " << data.SE_1["V2"] << " at Q = " << data.Q << endl;
  
  
  // need to iterate to determine MS bar mass for MChi to match equation (9) of Ibe et al.
  data.MChi = iterative_ms_bar_mass(data, "F7");
  
  data.P = data.MChi;
  data.do_tsil_all = true;
  //cout << "mw , mz = " << data.mw << ", " << data.mz << endl;
}
*/


// determine MSbar parameters

void MDM_spectrum::compute_spectra_MB_2loop()
{
	cout << "------------- Computing MB spectrum -------------" << endl;
  Self_energy se;
  long double tol = 1e-8;
 
  // determine MS bar input parameters at Q
  // only recompute the 1-loop functions for efficiency
  
  data.do_tsil_all = false;
  //data.mw = 79.06859081; // match the value FS changes this to
  long double mw_pole = data.mw;
  long double mz_pole = data.mz;
  long double MChi = data.MChi;
  long double alpha_SM = data.alpha;
  double diff = 1000;
  
  long double mwMS = 0.0;
  long double mzMS = 0.0;
  long double MChiMS = 0.0;
  long double alphaMS = 0.0;
  
  data.mz = 88.47644326;
  data.mw = 76.45989508; 
	
	diff = 100;
	double mwMS_previous = 0;
  while (diff > tol)
  {
		data.P = mz_pole;
		se.run_tsil(data);
		mzMS = pow( pow(mz_pole,2) - real(data.SE_1["V2"]),0.5 );

		cout << "SE( Z ) = " << data.SE_1["V2"] << " mzMS = " << mzMS << endl;
		
		data.P = mw_pole;	
		se.run_tsil(data);
		mwMS  = pow( pow(mw_pole,2) - real(data.SE_1["V3"]),0.5 );
		
		cout << "SE( W ) = " << data.SE_1["V3"] << " mwMS = " << mwMS << endl;
		
		/*
		data.P = data.Q;
		alphaMS = alpha_SM*( 1.0-real(extra_TSIL_interface_MDM::gammagamma_Chi(data,data.Q)) /pow(data.Q,2) );
		
		data.alpha = alphaMS;
		*/
		
		diff = abs(mzMS - data.mz) + abs(mwMS - mwMS_previous);
		cout << "diff = " << diff << endl;
		
		//data.mw = mwMS;
		mwMS_previous = mwMS;
		data.mz = mzMS;
	}
	
			// matching of SM to Wino model
	diff = 100;
	data.mw = mwMS_previous;
	
	
	data.P = MChi;
	
	data.do_tsil_all = true;
	
	
	
	
	std::cout << "MB: MZ pole = " << mz_pole << " ms_bar = " << data.mz << std::endl;
	std::cout << "MB: MW pole = " << mw_pole << " ms_bar = " << data.mw << std::endl;
  
  cout << "MB: alpha = " << data.alpha << endl;
	
  
}





bool MDM_spectrum::compute_spectra_flexiblesusy()
{
	
	cout << "------------- Computing FS spectrum -------------" << endl;
	
	Spectrum_generator_settings spectrum_generator_settings;
  
  QedQcd oneset;
  MDM_input_parameters input;
  MDM_spectrum_generator<Two_scale> spectrum_generator;
  
  MDM_slha_io slha_io;
  
  const std::string slha_file="flexiblesusy/models/MDM/LesHouches.in.MDM";
  
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
	const Problems<MDM_info::NUMBER_OF_PARTICLES>& problems
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
	MDM_slha<Two_scale> model(spectrum_generator.get_model());
	  
  // update data struct with computed spectrum
  
  // MS bar masses
  
  double g1 = pow(3./5.,0.5)*model.get_g1();
  double g2 = model.get_g2();
  
  // need to set masses to pole masses for this calculation along with the correct
  // vev and gauge couplings (forcing tree-level relationships)

  
  
  if (data.do_tsil_all)
  {
		
	  data.mw = model.get_MVWp();
	  data.mz = model.get_MVZ();
	  data.mh = model.get_Mhh();
	  data.mu = model.get_MFu(0);
	  data.mc = model.get_MFu(1);
	  //data.mt = model.get_MFu(2);
	  
	  data.md =  model.get_MFd(0);
	  data.ms =  model.get_MFd(1);
	  data.mb =  model.get_MFd(2);
	  
	  cout << "  u mass = " << model.get_MFu(0) << endl;
	  cout << "  c mass = " << model.get_MFu(1) << endl;
	  cout << "  t mass = " << model.get_MFu(2) << endl;
	  
	  cout << "  d mass = " << model.get_MFd(0) << endl;
	  cout << "  s mass = " << model.get_MFd(1) << endl;
	  cout << "  b mass = " << model.get_MFd(2) << endl;	  
	  
	  cout << "  e mass = " << model.get_MFe(0) << endl;
	  cout << "  mu mass = " << model.get_MFe(1) << endl;
	  cout << "  tau mass = " << model.get_MFe(2) << endl;	  
	  
	  
	  // check tree-level relationships:
	  
	  cout << "cw = mw/mz = " <<   data.mw/data.mz << endl;
	  
	  
	  
	  data.alpha = pow(g1*g2,2) / (4 * Pi * (g1*g1 + g2*g2));
	  cout << "alpha = " << data.alpha << endl;
	  
	  
	  //data.v = 2*data.mw/g2;
	  //data.v = model.get_v();
	  
		//double thetaW = ArcCos(Abs(model.get_ZZ(0,0)));
		//data.sw = Sin(thetaW);
		//data.cw = Cos(thetaW);
	  
	  //data.g1 = g1;
	  //data.g2 = g2;
	  
	  
	  //cout << "mw self energy FS = " << Re(model.self_energy_VWp(data.mw)) << ", MB = " << -real(data.SE_1["V3"]) << std::endl;
	}
	
	
	
	
  //data.g1 = g1;
  //data.g2 = g2;
  data.alpha = pow(g1*g2,2) / (4 * Pi * (g1*g1 + g2*g2));


 	data.SE_1["F7"] = model.get_MFn_pole_slha() - data.MChi;
  data.SE_1["F5"] = model.get_MFc_pole_slha() - data.MChi;
 
	
  return error;
}


void MDM_spectrum::compute_tsil()
{
	Self_energy se;
	cout << "run tsil " << endl;
	se.run_tsil(data);	
	
	// add derivatives of 1-loop self energies
	//data.SE_2["F7"] = data.SE_2["F7"] +  extra_TSIL_interface_MDM::F7_der(data);
	//data.SE_2["F6"] = data.SE_2["F6"] +  extra_TSIL_interface_MDM::F6_der(data);	
	//data.SE_2["F5"] = data.SE_2["F5"] +  extra_TSIL_interface_MDM::F5_der(data);	


	if (!data.do_tsil_all)
	{
		data.SE_2["F7"] = 0;
		data.SE_2["F6"] = 0;
		data.SE_2["F5"] = 0;
	}	
	
}


bool MDM_spectrum::compute_tsil_iterative()
{	
	Self_energy se;
	se.run_tsil(data);
	
	double m_F7 = iterative_ms_bar_mass(data,"F7");
	double m_F6 = 0.0L;
	
	if (m_F7 == 0)
	{
		return false;
	}
	
	if (m_F7 != 0)
	{
		m_F6 = iterative_ms_bar_mass(data,"F6_g1");
		if (m_F6 == 0)
		{
			return false;
		}
	}
	
	if ( (m_F7 == 0 ) || (m_F6 == 0 ) )
	{
		data.SE_1["F7"] = 0;
		data.SE_1["F6"] = 0;
		data.SE_2["F7"] = 0;
		data.SE_2["F6"] = 0;
	}
	else if (data.do_tsil_all)
	{
		data.SE_2["F7"] = m_F7 - data.MChi;
		data.SE_2["F6"] = m_F6 - data.MChi;
	}
	else
	{
		data.SE_1["F7"] = m_F7 - data.MChi;
		data.SE_1["F6"] = m_F6 - data.MChi;
	}	
return true;
}


double MDM_spectrum::get_deltam()
{
	return data.SE_1["F5"] - data.SE_1["F7"];
}

double MDM_spectrum::get_deltam_2loop()
{
	return data.SE_2["F5"] - data.SE_2["F7"];
}

double MDM_spectrum::get_deltam2()
{
	return data.SE_1["F6"] - data.SE_1["F7"];
}

double MDM_spectrum::get_deltam2_2loop()
{
	return data.SE_2["F6"] - data.SE_2["F7"];
}

double MDM_spectrum::get_double_charged_mass()
{
	return data.MChi + data.SE_1["F6_g1"] + data.SE_1["F6_g1"];
}

double MDM_spectrum::get_charged_mass()
{
	return data.MChi + data.SE_1["F5_g1"];
}

double MDM_spectrum::get_neutral_mass()
{
	return data.MChi + data.SE_1["F7_g1"] ;
}

double MDM_spectrum::get_charged_mass_2loop()
{
	return data.MChi + data.SE_1["F5_g1"]+data.SE_2["F5_g1"];
}

double MDM_spectrum::get_neutral_mass_2loop()
{
	return data.MChi + data.SE_1["F7_g1"] + data.SE_2["F7_g1"] ;
}
