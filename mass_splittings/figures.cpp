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
  int pts = 50;
  double M = 1000;
  int status = 0;
  
  double max_Q = 400; // (GeV)
  double min_Q = 50; // (GeV)
  data.MChi=M;
  data.P = M;
  
  double alpha_in = data.alpha;
  double MChi_in = data.MChi;
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  std::vector<double> mt={ 176.424,
    174.888,
    173.566,
    172.409,
    171.382,
    170.462,
    169.629,
    168.868,
    168.17,
    167.525,
    166.926,
    166.367,
    165.844,
    165.353,
    164.89,
    164.452,
    164.037,
    163.643,
    163.268,
    162.91,
    162.568,
    162.24,
    161.926,
    161.625,
    161.335,
    161.055,
    160.786,
    160.526,
    160.275,
    160.033,
    159.798,
    159.57,
    159.349,
    159.135,
    158.927,
    158.725,
    158.528,
    158.337,
    158.15,
    157.969,
    157.792,
    157.619,
    157.451,
    157.286,
    157.125,
    156.968,
    156.814,
    156.664,
    156.517,
    156.372,
    156.231};
  
  
  for (int i = 0; i < pts+1 ; i++)
  {
    
    MSSM mssm;
    
    data.mt = mt[i];
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    
    
    mssm.set_SM_parameters_1loop(data);
    data.P = data.MChi;
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    data.do_tsil_all = false;
    Self_energy se1;
    se1.run_tsil(data);
    data.do_tsil_all = true;
    
    double delta_m_1 = data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
    
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.P = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    mssm.set_SM_parameters_2loop(data);
    data.P = data.MChi;
    data.Q = pow((i)*(max_Q - min_Q)/pts + min_Q,2);
    
    Self_energy se2;
    se2.run_tsil(data);
    
    
    double delta_m_2 = (data.SE_1["F12_g1"] - data.SE_1["F11_g1"]) + ( data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - mssm.add_1loop_der(data) );
    
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
  double M = 0;
  int status = 0;
  
  double max_M = 4000; // (GeV)
  double min_M = 90; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
      
  double alpha_in = data.alpha;
  
  double mw_in = data.mw;
  double mz_in = data.mz;
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    
    M = pow(10,n);
    double MChi_in = M;
    data.MChi=M;
    data.P = M;

    
    
    data.mt = 163.3;
    
    data.Q = pow(163.3,2);
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    
    
    MSSM mssm;

    mssm.set_SM_parameters_1loop(data);
    data.P = data.MChi;
    data.Q = pow(163.3,2);
    
    data.do_tsil_all = false;
    Self_energy se1;
    se1.run_tsil(data);
    data.do_tsil_all = true;
    
    double delta_m_1 = data.SE_1["F12_g1"] - data.SE_1["F11_g1"];
    
    data.alpha = alpha_in;
    data.MChi = MChi_in;
    data.P = MChi_in;
    data.mw = mw_in;
    data.mz = mz_in;
    mssm.set_SM_parameters_2loop(data);
    data.P = data.MChi;
    data.Q = pow(163.3,2);
    
    Self_energy se2;
    se2.run_tsil(data);
    
    
    double delta_m_2 = (data.SE_1["F12_g1"] - data.SE_1["F11_g1"]) + ( data.SE_2["F12_g1"] - data.SE_2["F11_g1"] - mssm.add_1loop_der(data) );
    
    myfile << M << " " << delta_m_1 <<  " " << delta_m_2 << endl;
    
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
