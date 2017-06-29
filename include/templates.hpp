#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <unistd.h>


using namespace std;

namespace templates
{
  
  void self_energy_hpp_header(ofstream &file);
  
  void data_input_block(ofstream &file);
  
  void self_energy_src_preamble(ofstream &file, std::string particle_name, int subgroup);
  
  void amp_preamble(ofstream &file);
  
  void status_update(ofstream &file);
  
  void user_input_guide();
  
  void print_math_header(std::string &input);
  
}

#endif