#ifndef WRITE_TSIL_INI_H
#define WRITE_TSIL_INI_H

/* 

This class adds the ability to significantly reduce the number TSIL_Evaluate calls required during
a self energy calculation by taking advantage of the symmetries between basis integrals in an entirely
generic way

Sep 2016

*/


#include "bases.hpp"
#include "utils.hpp"

using namespace std;
using namespace utils;

class Print_dotsil
{

  vector<string> names;
  std::map<std::string, Bases> base_map;

public:
  
  Print_dotsil(vector<string> names, std::map<std::string, Bases> base_map) : names(names) , base_map(base_map) {}
  
  void print_to_file(ofstream &myfile);
  void sort_integrals();
  
  void Print_dotsil::get_poss_eval(Bases base);
  
  
};

#endif