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


#ifndef PI
#define PI 4.0L*atan(1.0L)
#endif

using namespace std;

namespace extra_TSIL_interface
{
#include TSIL_PATH

  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}
  double d1loop(Data &data)
  {
    double Pi = PI;
    double g = data.g;
    double C = g*g/(16.*Pi*Pi);
    double Ms2 = TSIL_POW(data.Ms,2);
    
    TSIL_COMPLEXCPP result = -0.5L*C*TSIL_dBds_(Ms2,Ms2,TSIL_POW(data.P,2),data.Q);
    
    return real(result);
  }
}


double iterative_mass(Data data, int loop_order)
{
  double M_tree = data.Ms;
  double M_pole = data.Ms;
  data.P = M_tree;
  double diff = 1;
  double precision = 1e-30;
  int iteration =0;
  do
  {
    Self_energy se;
    se.run_tsil(data);
    
    if (loop_order == 1)
    {
      M_pole = TSIL_POW(M_tree*M_tree + (data.SE_1["S1"]),0.5);
    }
    
    if (loop_order == 2)
    {
      M_pole = TSIL_POW(M_tree*M_tree + (data.SE_1["S1"]+data.SE_2["S1"]),0.5);
    }
    
    diff = abs(M_pole - data.P);
    
    data.P = M_pole;
    
    iteration++;
    //cout << "diff = " << diff << endl;
    
  } while (diff > precision  && iteration < 500);
  
  if (iteration == 500)
  {
    cout << "pole mass did not converge" << endl;
  }
  
  return M_pole;
}



int main()
{
  // input data (see other examples for more advanved input and user interface)
  Data data;
  data.P = 10.;
  data.Q = 100.;
  data.Ms = 10.;
  data.lambda = 0.01;
  data.g = 0.5;
  
  // initialise self energy class
  Self_energy se;
  se.run_tsil(data);
  
  // request self energies for particle "S1"
  double one_loop = data.SE_1["S1"];
  double two_loop = one_loop*(1.+extra_TSIL_interface::d1loop(data)) + data.SE_2["S1"];
  double two_loop_noderivative = one_loop + data.SE_2["S1"];
  
  
  
  cout.precision(17);
  
  cout << "--- Explicit pole masses --- " << endl;
  cout << "One-loop pole mass = " << pow(data.Ms*data.Ms + one_loop,0.5) << endl;
  cout << "Two-loop pole mass = " << pow(data.Ms*data.Ms + two_loop,0.5) << endl;
  
  cout << "--- Explicit pole masses without derivative term --- " << endl;
  cout << "Two-loop pole mass = " << pow(data.Ms*data.Ms + two_loop_noderivative,0.5) << endl;
  
  cout << "--- Iterative pole masses --- " << endl;
  
  cout << "One-loop pole mass = " << iterative_mass(data,1) << endl;
  cout << "Two-loop pole mass = " << iterative_mass(data,2) << endl;
  
  //cout << "derivative term = " << one_loop*extra_TSIL_interface::d1loop(data) << endl;
  
  cout << "--- Differences --- " << endl;
  
  cout << "One-loop " << endl;
  
  cout << "Explicit - Implicit = " << pow(data.Ms*data.Ms + one_loop,0.5) -  iterative_mass(data,1) << endl;
  
  cout << "Two-loop " << endl;
  
  cout << "Explicit - Implicit = " << pow(data.Ms*data.Ms + two_loop,0.5) -  iterative_mass(data,2) << endl;
  
  cout << "Explicit no derivative - Implicit = " << pow(data.Ms*data.Ms + two_loop_noderivative,0.5) -  iterative_mass(data,2) << endl;
  
  cout << "Explicit no derivative - Explicit = " << pow(data.Ms*data.Ms + two_loop_noderivative,0.5) -  pow(data.Ms*data.Ms + two_loop,0.5) << endl;
  
  
  return 0;
}


















