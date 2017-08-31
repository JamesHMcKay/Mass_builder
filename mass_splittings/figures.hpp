#ifndef FIGURES_H
#define FIGURES_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "data.hpp"

#include "self_energy.hpp"
#include "cmake_variables.hpp"

#include "MSSM.hpp"
#include "decays.hpp"

class Figures
{
private:
  
public:
  
  Figures(){};
  
	void plot_Q(Data data);
	
	void plot_M(Data data);
	
	void plot_M_flexiblesusy(Data data);
	
	void plot_uncertainties(Data data);
	
	void test(Data data);
	
	void plot_M_flexiblesusy_2loop(Data data);
  
};

#endif
