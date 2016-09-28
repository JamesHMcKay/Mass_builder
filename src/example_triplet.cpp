/*
 Mass Builder 
 
 James McKay
 Sep 2016
 
 -- example_triplet.cpp --
 
 compute full two-loop self energy including derivative of 1-loop functions
 
 requires an input list flag at runtime: ./triplet -i models/MDM/input.txt
 */

#include "data.hpp"
#include "calc_amplitudes.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"
#include "supplements.hpp"

using namespace std;
using namespace supplementary_code;
using namespace utils;

Self_energy se;

double pole_mass_F5(Data data)
{
  double Mp = data.MChi + (data.SE_1["F5"]+data.SE_2["F5"]);
  return Mp;
}

double pole_mass_F6(Data data)
{
  double Mp = data.MChi + (data.SE_1["F6"]+data.SE_2["F6"]);
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

  
  se.run_tsil(data);
  
  double se_two_loop = data.SE_2["F5"];
  cout << "self energy from two loop only = " << se_two_loop << endl;
  Supplements supp(data);
  supp.add_derivatives(data);
  cout << "derivative terms = " << data.SE_2["F5"] - se_two_loop << endl;
  
  
  
  
  
  return 0;
}
