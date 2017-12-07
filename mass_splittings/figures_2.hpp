#ifndef FIGURES_2_H
#define FIGURES_2_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "data.hpp"

#include "self_energy.hpp"
#include "cmake_variables.hpp"

#include "ew_triplet_spectrum.hpp"
#include "mdm_spectrum.hpp"


#include "decays.hpp"
#include "utils.hpp"


struct Limits
{
	double upper;
	double lower;
	
	double upper_Q;
	double lower_Q;
	
};

template <class T>
class Figures_2
{
private:
  
public:
  
  Figures_2(){};
  
	void two_loop_plots(Data data, string group);
	
	void two_loop_plots_uncertainties(Data data, string group);
	
	void test_MDM(Data data);
	
	//Limits find_uncertainties(Data data);
	
	//Limits find_uncertainties2(Data data);
	
};



#endif
