#ifndef FIGURES_H
#define FIGURES_H

#include "pv.hpp"
#include "data.hpp"
#include "decays.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>

using namespace std;
typedef complex<long double> dcomp;

class Figures
{
  private:
  
  Data data;

  public:
  
  Figures(Data _data){data=_data;}
  
  void plot_mass_splittings();
  void plot_pole_masses_n();
  void plot_pole_masses_c();
  void plot_SigmaK();
  void plot_limits();
  void plot_approx();
  void plot_derivatives();
  long double theory(long double M);
  
  void plot_decay();
  void plot_mass_unc();
  void plot_mass_splitting_unc();
  
};

#endif