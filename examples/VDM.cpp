/*
 Mass Builder 
 
 James McKay
 May 2017
 
 -- VDM.cpp --
 
 determine mass splittings for vector dark matter triplet model
 
 requires an input list flag at runtime: ./VDM -i models/VDM/input.txt
 */

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"
  
#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

double Pi;

using namespace std;


namespace extra_TSIL_interface
{
  #include TSIL_PATH
  using namespace std;
  TSIL_DATA bar;
 
#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif
// define subroutines here
  TSIL_REAL Q2;
  TSIL_REAL p;
  bool exclude_photon_pole = false;
  TSIL_COMPLEXCPP Log(TSIL_REAL a){
  complex<double> s(a/Q2,-0.000);return log(s);}
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
 
  TSIL_COMPLEXCPP Ah, Ap, Av, Aw, Az ;
  TSIL_COMPLEXCPP Bnp, Bph, Bvh, Bwp, Bwv ;
  TSIL_COMPLEXCPP Bzp, Bpn, Bhp, Bhv, Bpw ;
  TSIL_COMPLEXCPP Bvw, Bpz ;
  TSIL_COMPLEXCPP  i;
  
  double Zeta;
  TSIL_REAL  null,  null2 ,  MZ,  MZ2 ,  MW,  MW2 ,  MVp,  MVp2 ,  MV0,  MV02 ,  MH, MH2 ;

  TSIL_REAL Pi;
  TSIL_REAL CW, CW2, EE, MassBuilderCTM1, MassBuilderCTM2, MassBuilderCTZ1, MassBuilderCTZ2, SW, SW2, a ;

void DoTSIL(Data data)
{
  SW = data.SW,   CW2 = data.CW2,   CW = data.CW,   SW2 = data.SW2,   EE = data.EE,   a = data.a ;

   TSIL_REAL Q2 = pow(data.Q,2);
   TSIL_REAL s = pow(data.P,2);
   p = data.P;
   MassBuilderCTM1 = 0;
   MassBuilderCTZ1 = 0;
   MassBuilderCTM2 = 0;
   MassBuilderCTZ1 = 0;
  null = data.null, null2 = TSIL_POW(data.null, 2) ,   MZ = data.MZ, MZ2 = TSIL_POW(data.MZ, 2) ,   MW = data.MW, MW2 = TSIL_POW(data.MW, 2) ,   MVp = data.MVp, MVp2 = TSIL_POW(data.MVp, 2) ,   MV0 = data.MV0, MV02 = TSIL_POW(data.MV0, 2) ,   MH = data.MH , MH2 = TSIL_POW(data.MH, 2) ;

    dcomp ii=-1;ii=sqrt(ii);i=ii;
   Pi=PI;
    SW =   TSIL_POW(SW2,0.5);
    CW2 =  1.0-SW2;
    CW =  TSIL_POW(CW2,0.5);
    
    TSIL_REAL a = 1., b = 1.;
    
  Ah = -a*TSIL_A_ (MH2 , Q2);

  Ap = -a*TSIL_A_ (MVp2 , Q2);

  Av = -a*TSIL_A_ (MV02 , Q2);

  Aw = -a*TSIL_A_ (MW2 , Q2);

  Az = -a*TSIL_A_ (MZ2 , Q2);

  Bnp = b*TSIL_B_ (null2, MVp2, s, Q2);
  Bpn = Bnp;

  Bph = b*TSIL_B_ (MVp2, MH2, s, Q2);
  Bhp = Bph;

  Bvh = b*TSIL_B_ (MV02, MH2, s, Q2);
  Bhv = Bvh;

  Bwp = b*TSIL_B_ (MW2, MVp2, s, Q2);
  Bpw = Bwp;

  Bwv = b*TSIL_B_ (MW2, MV02, s, Q2);
  Bvw = Bwv;

  Bzp = b*TSIL_B_ (MZ2, MVp2, s, Q2);
  Bpz = Bzp;
  
  
  
  
  }
  
  double one_loop_analytic(Data data)
  {
		DoTSIL(data);
		
		Pi = PI;
    p = data.P;
    EE = data.EE;
    TSIL_REAL SW2 = data.SW2;
    SW =   TSIL_POW(SW2,0.5);
    
    TSIL_REAL a = data.a;
		
		TSIL_REAL MChi = data.MVp;
		
		
		TSIL_COMPLEXCPP SEc = (-30*Ap*Power(EE,4)*Power(MChi,2) + 15*Aw*Power(EE,4)*Power(MChi,2) + 15*Az*Power(EE,4)*Power(MChi,2) + 
     15*Bwp*Power(EE,4)*Power(MChi,4) + 15*Bzp*Power(EE,4)*Power(MChi,4) + 15*Bwp*Power(EE,4)*Power(MW,4) + 
     15*Bzp*Power(EE,4)*Power(MZ,4) + 42*Ap*Power(EE,4)*Power(p,2) + 21*Aw*Power(EE,4)*Power(p,2) + 21*Az*Power(EE,4)*Power(p,2) - 
     60*Power(EE,4)*Power(MChi,2)*Power(p,2) - 48*Bwp*Power(EE,4)*Power(MChi,2)*Power(p,2) - 
     48*Bzp*Power(EE,4)*Power(MChi,2)*Power(p,2) - 4*Power(EE,4)*Power(p,4) - 57*Bwp*Power(EE,4)*Power(p,4) - 
     57*Bzp*Power(EE,4)*Power(p,4) + 3*Power(EE,4)*Power(MW,2)*
      (5*Ap - 5*Aw - 2*(5*Bwp*Power(MChi,2) + 5*Power(p,2) + 8*Bwp*Power(p,2))) + 
     3*Power(EE,4)*Power(MZ,2)*(5*Ap - 5*Az - 2*(5*Bzp*Power(MChi,2) + 5*Power(p,2) + 8*Bzp*Power(p,2))) - 
     15*Bzp*Power(EE,4)*Power(MZ,4)*Power(SW,2) + 3*Power(EE,4)*Power(MZ,2)*
      (-5*Ap + 5*Az + 2*(5*Bzp*Power(MChi,2) + 5*Power(p,2) + 8*Bzp*Power(p,2)))*Power(SW,2) - 
     3*Power(EE,2)*(5*Bzp*Power(EE,2)*Power(MChi,4) - 6*a*Ah*Power(p,2) - 12*a*Aw*Power(p,2) - 
        16*Bzp*Power(EE,2)*Power(MChi,2)*Power(p,2) - 19*Bzp*Power(EE,2)*Power(p,4) + 
        Bnp*Power(EE,2)*(-5*Power(MChi,4) + 16*Power(MChi,2)*Power(p,2) + 19*Power(p,4)) + 
        Az*(-6*a*Power(p,2) + Power(EE,2)*(5*Power(MChi,2) + 7*Power(p,2))))*Power(SW,2) + 
     288*Power(a,2)*Bhp*Power(MW,2)*Power(p,2)*Power(SW,4))/(288.*Power(EE,2)*Power(p,2)*Power(Pi,2)*Power(SW,2));
     
     TSIL_COMPLEXCPP SEn = (15*Aw*Power(EE,4)*Power(MChi,2) + 15*Bwp*Power(EE,4)*Power(MChi,4) + 15*Bwp*Power(EE,4)*Power(MW,4) + 21*Aw*Power(EE,4)*Power(p,2) - 
     30*Power(EE,4)*Power(MChi,2)*Power(p,2) - 48*Bwp*Power(EE,4)*Power(MChi,2)*Power(p,2) - 2*Power(EE,4)*Power(p,4) - 
     57*Bwp*Power(EE,4)*Power(p,4) - 3*Ap*Power(EE,4)*(5*Power(MChi,2) - 5*Power(MW,2) - 7*Power(p,2)) - 
     3*Power(EE,4)*Power(MW,2)*(5*Aw + 2*(5*Bwp*Power(MChi,2) + 5*Power(p,2) + 8*Bwp*Power(p,2))) + 
     9*a*(Ah + 2*Aw + Az)*Power(EE,2)*Power(p,2)*Power(SW,2) + 144*Power(a,2)*Bhp*Power(MW,2)*Power(p,2)*Power(SW,4))/
   (144.*Power(EE,2)*Power(p,2)*Power(Pi,2)*Power(SW,2));
		
		
		TSIL_COMPLEXCPP result1 = 0.5L*(SEc-SEn)/MChi;
		/*
		TSIL_COMPLEXCPP result2 = -(1.L/8.L)*(  pow(real(SEc),2) - pow(real(SEn),2) ) /pow(MChi,3);
		
		TSIL_COMPLEXCPP result3 = -(1.L/16.L)*(  pow(real(SEc),3) - pow(real(SEn),3) ) /pow(MChi,5);
		
		TSIL_COMPLEXCPP result4 = -(5.L/128.L)*(  pow(real(SEc),4) - pow(real(SEn),4) ) /pow(MChi,7);

		// compute difference between pole mass and expansion
		
		double Mpole_expansion = MChi + 0.5L*real(SEn)/MChi ;
		
		double Mpole = pow( pow(MChi,2) + real(SEn) , 0.5 ) ;
		
		cout <<  "M = " << MChi << ", Mpole = " << Mpole << " , " << Mpole_expansion << " , diff = " << abs(Mpole_expansion-Mpole)  <<  ", SEn/M^2 = " << real(SEn/pow(MChi,2)) << endl;
		
		
		//cout << "M = " << MChi <<  ", SEn/M^2 = " << real(SEn/pow(MChi,2)) << ", SEn^2/M^4 " << pow(real(SEn),2) / pow(MChi,4) << endl;
		//cout << "M = " << MChi <<  ", SEc/M^2 = " << real(SEc/pow(MChi,2)) << ", SEc^2/M^4 " << pow(real(SEc),2) / pow(MChi,4) << endl;
		*/
		
		return real(-result1);
	}
}



double pole_mass_V5(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  double Mp = pow(data.MVp*data.MVp + (data.SE_1["V5"]),0.5);
  
  if (data.MVp*data.MVp + (data.SE_1["V5"]) < 0)
  {
    cout << "tachyon mass!" << endl;
  }
  
  
  return Mp;
}

double pole_mass_V6(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  double Mp = pow(data.MV0*data.MV0 + (data.SE_1["V6"]),0.5);
  
  if (data.MVp*data.MVp + (data.SE_1["V6"]) < 0)
  {
    cout << "tachyon mass!" << endl;
  }
  return Mp;
}

double pole_mass_V5_expansion(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  
  double Mp = data.MV0 + 0.5L*data.SE_1["V5"] / data.MV0;

  return Mp;
}



double pole_mass_V6_expansion(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  
  double Mp = data.MV0 + 0.5L*data.SE_1["V6"] / data.MV0;

  return Mp;
}



// determine MSbar parameters
void compute_spectra(Data &data)
{
  double alpha = pow(data.EE,2) / (4.*Pi);
  
  //cout << "alpha (mz) = " << alpha ;
  double A = (7./(30.*Pi)); // SM
  double alpha_mz = alpha;
  double mu0 = data.MZ;
  double mu = data.Q;
  
  alpha = pow(   1.0/alpha_mz -  A * log( mu/mu0) , -1);
  //cout << ", alpha(" << data.Q << ") = " << alpha << endl;
  
  data.EE = pow( (4.*Pi) * alpha , 0.5) ;
}

int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy se;
  Data data(options);
  ofstream deltam;
  deltam.open ("models/VDM/output/mass_splittings.txt");
  
    // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(10.0L);
  
  std::vector<double> Q(5);
  double EE_mz = data.EE;
  Pi = PI;
  // number of points to plot
  int pts = 100;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    data.MVp = M;
    data.MV0 = M;
    data.P = M;
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    deltam << M ;
    
    for (int i = 0; i < 5 ; i++)
    {
			data.Q = Q[i];
			compute_spectra(data);
	    deltam << " " << pole_mass_V5(data) - pole_mass_V6(data);
	    data.EE = EE_mz;
		}
		
		data.Q = 173.15;
		deltam << " " << extra_TSIL_interface::one_loop_analytic(data);
		
		//deltam << " " << pole_mass_V5_expansion(data) - pole_mass_V6_expansion(data);
		
		deltam << endl;
  }
  
  
  cout << "VDM mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_VDM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  
	deltam.close();
  
  return 0;
}
