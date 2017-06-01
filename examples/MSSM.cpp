/*
 Mass Builder
 
 -- MSSM.cpp --
 
 example application, determing mass splitting between neutral and charged winos in MSSM
 
 requires an input list flag at runtime: ./MSSM -i models/MSSM/input.txt
 */

#include "data.hpp"
#include "calc_amplitudes.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"

using namespace std;
using namespace utils;


double iterative_mass_F5(Data data)
{
  
  double M = data.MChi;
  double Mp = data.MChi;
  
  double MFn=M,M_tree=M,new_MFn,old_MFn=M,p;
  
  double diff = 1;
  double precision = 0.0001;
  int iteration =0;
  
  do{
    p=old_MFn;
    
    data.P = p;
    Self_energy se;
    se.run_tsil(data);
    
    double M_1loop=M_tree + (data.SE_1["F12_g1"]+data.SE_2["F12_g1"]);
    MFn=M_1loop;
    new_MFn=MFn;
    diff=abs(new_MFn-old_MFn);
    old_MFn=new_MFn;
    iteration++;
    //cout << "diff = " << diff << endl;
    
  } while (diff > precision  && iteration < 500);
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  Mp = new_MFn;
  
  //cout << "----- done ----- " << endl;
  return Mp;
}

double iterative_mass_F6(Data data)
{
  
  double M = data.MChi;
  double Mp = data.MChi;
  
  double MFn=M,M_tree=M,new_MFn,old_MFn=M,p;
  
  double diff = 1;
  double precision = 0.0001;
  int iteration =0;
  
  do{
    p=old_MFn;
    
    data.P = p;
    Self_energy se;
    se.run_tsil(data);
    
    double M_1loop=M_tree +(data.SE_1["F11_g1"]+data.SE_2["F11_g1"]);
    MFn=M_1loop;
    new_MFn=MFn;
    diff=abs(new_MFn-old_MFn);
    old_MFn=new_MFn;
    iteration++;
    //cout << "diff = " << diff << endl;
    
  } while (diff > precision && iteration < 500);
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  //cout << "----- done ----- " << endl;
  
  Mp = new_MFn;
  return Mp;
}


double pole_mass_F5(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  //Supplements supp(data);
  //supp.add_derivatives(data);
  double Mp = data.MChi + (data.SE_1["F12_g1"]+data.SE_2["F12_g1"]);
  return Mp;
}

double pole_mass_F6(Data data)
{
  Self_energy se;
  se.run_tsil(data);
  //Supplements supp(data);
  //supp.add_derivatives(data);
  double Mp = data.MChi + (data.SE_1["F11_g1"]+data.SE_2["F11_g1"]);
  return Mp;
}



int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy se;
  Data data(options);
  
  ofstream myfile;
  myfile.open ("models/MSSM/output/mass_splittings.txt");
  int pts = 30;
  double n = 0;
  double M=0;
  int status = 0;
  for (int i = 0; i < pts ; i++)
  {
    n=(float(i)/float(pts))*5;
    M= pow(10,n);
    data.MChi=M;
    data.P = M;
    
    double delta_m_it=iterative_mass_F5(data) - iterative_mass_F6(data);
    data.MChi=M;
    data.P = M;
    M= pow(10,n);
    double delta_m=pole_mass_F5(data) - pole_mass_F6(data);
    
    myfile << M << " " << delta_m_it << " " << delta_m << endl;
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
  
  return 0;
}
