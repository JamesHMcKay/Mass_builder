/*
  Mass Builder

  James McKay
  May 2017

  -- Scalar.cpp --

  example routine to access self energies
*/

#include "data.hpp"
#include "compute_amp.hpp"
#include "self_energy.hpp"

using namespace std;

int main()
{
  // input data (see other examples for more advanved input and user interface)
  Data data;
  data.P = 10.;
  data.Q = 100.;
  data.Ms = 10.;
  data.lambda = 0.1;
  data.g = 0.1;
  
  // initialise self energy class
  Self_energy se;
  se.run_tsil(data);
  
  // request self energies for particle "S1"
  double one_loop = data.SE_1["S1"];
  double two_loop = data.SE_2["S1"];
  
  cout << "One-loop self energy = " << one_loop << endl;
  cout << "Two-loop self energy = " << two_loop << endl;
  cout << "Total self energy = " << one_loop+two_loop << endl;
  
  return 0;
}