/*
 Mass Builder
 
 -- splittings.cpp --

 Master routine for calling plotting functions for mass
 splitting results
 
 */

#include "data.hpp"
#include "self_energy.hpp"
#include "cmake_variables.hpp"
#include "figures.hpp"
#include "MSSM.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  User_input user(argc,argv);
  user.user_interface();
  Options options = user.options;
  
  if (options.input_list == "") {cout << "please enter an input list" << endl; return 0;}
  
  Data data(options);
  
  Figures fig;
  
  fig.plot_M(data);
  
  return 0;
}
