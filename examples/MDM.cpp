/*
 Mass Builder
 
 -- MSSM.cpp --
 
 example application, determing mass splitting between neutral and charged winos in MSSM
 
 requires an input list flag at runtime: ./MSSM -i models/MSSM/input.txt
 */


#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"

using namespace std;

#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif


double iterative_mass(Data data, string particle,int loop_order)
{
  double M_tree = data.MChi;
  double M_pole = data.MChi;
  data.P = M_tree;
  double diff = 1;
  double precision = 1e-4;
  int iteration =0;
  int loop_order_temp = loop_order;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    if (iteration==0)
    {
     loop_order = 1;
    }
    
    if (loop_order == 1)
    {
      M_pole = M_tree + (data.SE_1[particle]);
    }
    
    if (loop_order == 2)
    {
      M_pole = M_tree + (data.SE_1[particle]+data.SE_2[particle]);
    }
    
    diff = abs(M_pole - data.P);
    
    data.P = M_pole;
    
    if (iteration == 0)
    {
      loop_order = loop_order_temp;
    }
    
    iteration++;
    cout<< "\r" << "M_pole - p = " << diff << " GeV";
    std::cout << std::flush;
  } while (diff > precision  && iteration < 500);
  
  cout<< "\r" << "M_pole - p = " << diff << " GeV";
  cout << "\n";
  
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  return M_pole;
}








void plot_M(Data data)
{
  
  
  ofstream myfile;
  myfile.open ("models/MDM/output/mass_splittings.txt");
  int pts = 30;
  double n = 0;
  double M = 0;
  int status = 0;
  
  double max_M = 10000; // (GeV)
  double min_M = 10; // (GeV)
  
  double logMax = log10(max_M);
  double logMin = log10(min_M);
  
  
  for (int i = 0; i < pts+1 ; i++)
  {
    n = (i)*(logMax - logMin)/pts + logMin;
    
    
    
    M = pow(10,n);
    data.MChi=M;
    data.P = M;
    
    data.Q = M;
    
    // compute explicit mass splitting
    
    Self_energy se;
    se.run_tsil(data);
    
    double delta_m1_explicit_1loop = (data.SE_1["F5"]) - (data.SE_1["F7"]);
    
    double delta_m2_explicit_1loop = (data.SE_1["F6"]) - (data.SE_1["F7"]);
    
    
    
    // compute iterative mass splitting
    
    
    
    double delta_m1_iterative_1loop = iterative_mass(data,"F5",1) - iterative_mass(data,"F7",1);
    
    double delta_m2_iterative_1loop = iterative_mass(data,"F6",1) - iterative_mass(data,"F7",1);
    
    
    myfile << M << " " << delta_m1_iterative_1loop <<  " " << delta_m2_iterative_1loop << " " << delta_m1_explicit_1loop <<  " " << delta_m2_explicit_1loop << endl;
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









int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  Data data(options);
  
  plot_M(data);
  
  return 0;
}
