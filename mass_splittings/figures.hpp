#ifndef FIGURES_H
#define FIGURES_H

#include "data.hpp"

#include "self_energy.hpp"
#include "cmake_variables.hpp"

#include "MSSM.hpp"

class Figures
{
private:
  
public:
  
  Figures(){};
  
	void plot_Q(Data data);
	
	void plot_M(Data data);
	
	void plot_M_flexiblesusy(Data data);
  
};

#endif
