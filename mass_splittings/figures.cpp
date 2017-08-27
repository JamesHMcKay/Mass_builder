/*
 Mass Builder
 
 -- figures.cpp --
 
 Contains functions for generating the data for plots
 
 */


#include "figures.hpp"

using namespace std;

void Figures::plot_Q(Data data)
{
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings_Q.txt");
  int pts = 15;
  double M = 1000;
  int status = 0;
  
  double max_Q = 400; // (GeV)
  double min_Q = 50; // (GeV)
  data.MChi = M;
  data.P = M;
  
  for (int i = 0; i < pts+1 ; i++)
  {
    
    data.mt = 163.3;  
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    MSSM mssm_1loop(data);
    mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    
    double delta_m_1 = mssm_1loop.data.SE_1["F12_g1"] - mssm_1loop.data.SE_1["F11_g1"];
    
    
    MSSM mssm_2loop(data);
    mssm_2loop.compute_spectra_MB_2loop();
    mssm_2loop.compute_tsil();
    
    double delta_m_2 = (mssm_2loop.data.SE_1["F12_g1"] - mssm_2loop.data.SE_1["F11_g1"]) + ( mssm_2loop.data.SE_2["F12_g1"] - mssm_2loop.data.SE_2["F11_g1"]);
    
    myfile << pow(data.Q,0.5) << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  cout << "mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_MSSM_Q.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  myfile.close();  
}

void Figures::plot_M(Data data)
{
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 10;
  double n = 0;
  int status = 0;
  
  double max_M = 4000; // (GeV)
  double min_M = 90; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
      
  
  data.mt = 163.3;  
  data.Q = pow(163.3,2);
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    
    data.MChi = pow(10,n);
    data.P = data.MChi;
    
    MSSM mssm_1loop(data);
    mssm_1loop.compute_spectra_MB_1loop();
    mssm_1loop.compute_tsil();
    
    double delta_m_1 = mssm_1loop.data.SE_1["F12_g1"] - mssm_1loop.data.SE_1["F11_g1"];
    
    MSSM mssm_2loop(data);
    mssm_2loop.compute_spectra_MB_2loop();
    mssm_2loop.compute_tsil();
    
    double delta_m_2 = (mssm_2loop.data.SE_1["F12_g1"] - mssm_2loop.data.SE_1["F11_g1"]) + ( mssm_2loop.data.SE_2["F12_g1"] - mssm_2loop.data.SE_2["F11_g1"]);
    
    myfile << data.MChi << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    
    status=(float(i)/pts)*100;
    cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  
  status=100;
  cout<< "\r" << "computing mass splittings . . . " << status << "% complete ";
  cout << "\n";
  
  cout << "example mass splitting routine complete" << endl;
  cout << "now run: "<< endl;
  cout << "          python examples/plot_MSSM.py "<< endl;
  cout << "to make plot in this directory "<< endl;
  
  myfile.close();
}



void Figures::plot_M_flexiblesusy(Data data)
{
  // plot MDM mass splittings
  
  ofstream myfile_a;
  myfile_a.open ("mass_splittings/mass_splittings_iterative.txt");
  
  //ofstream myfile_b;
  //myfile_b.open ("mass_splittings/mass_splittings_explicit.txt");
  
  // set range of plot
  long double logMax = log10(1.0e4L);
  long double logMin = log10(100.0L);
  
  std::vector<double> Q(5);
  
  // number of points to plot
  int pts = 30;
  for (int i=0;i<pts+1;i++)
  {
    double n = i*(logMax - logMin)/pts + logMin;
    double M = pow(10.0L,n);
    
    data.MChi = M;
    data.P = data.MChi;
    
    Q[0] = 2.0 * 173.15;
    Q[1] = 0.5 * 173.15;
    Q[2] = 0.5 * M ;
    Q[3] = M ;
    Q[4] = 2.0 * M;
    
    myfile_a << M ;
    
    for (int i = 2; i < 3 ; i++)
    {
			data.Q = pow(Q[i],2);
			MSSM mssm(data);
			//mssm.compute_spectra_flexiblesusy();
			mssm.compute_spectra_MB_2loop();
      mssm.compute_tsil_iterative();
			
			myfile_a << " " << mssm.data.SE_2["F12_g1"] - mssm.data.SE_2["F11_g1"];
		}
    
    myfile_a << endl;
    
  }
  myfile_a.close();
  //myfile_b.close();
  //system("python ../Figures/mass_splittings.py");
}
