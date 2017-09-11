#ifndef FIGURES_H
#define FIGURES_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "data.hpp"

#include "self_energy.hpp"
#include "cmake_variables.hpp"

#include "mssm_spectrum.hpp"

#include "ew_triplet_spectrum.hpp"
#include "mdm_spectrum.hpp"


#include "decays.hpp"

template <class T>
class Figures
{
private:
  
public:
  
  Figures(){};
  
	void plot_Q(Data data);
	
	void plot_M(Data data);
	
	void plot_M_flexiblesusy(Data data);
	
	void plot_uncertainties(Data data);
	
	void plot_M_2loop_iterative(Data data);
	
	void plot_M_2loop_explicit(Data data);
  
};



#endif
