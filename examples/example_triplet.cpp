/*
 Mass Builder 
 
 James McKay
 Sep 2016
 
 -- example_triplet.cpp --
 
 compute full two-loop self energy including derivative of 1-loop functions
 
 requires an input list flag at runtime: ./triplet -i models/EW_triplet/input.txt
 */

#include "data.hpp"
#include "calc_amplitudes.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"
#include "supplements.hpp"

using namespace std;
using namespace supplementary_code;
using namespace utils;

Self_energy self_energy;

double pole_mass_F1(Data data)
{
  double Mp = data.MChi + (data.SE_1["F1"]+data.SE_2["F1"]);
  return Mp;
}

double pole_mass_F2(Data data)
{
  double Mp = data.MChi + (data.SE_1["F2"]+data.SE_2["F2"]);
  return Mp;
}



int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Self_energy self_energy;
  Data data(options);

  self_energy.run_tsil(data);
  
  cout << "one-loop mass splitting = " << data.SE_1["F2"] - data.SE_1["F1"] << endl;

  cout << "two-loop mass splitting = " << data.SE_2["F2"] - data.SE_2["F1"] << endl;
  
  Supplements supplemets(data);
  supplemets.add_derivatives(data);

  cout << "two-loop mass splitting including derivative terms = " << data.SE_2["F2"] - data.SE_2["F1"] << endl;
  
  return 0;
}
