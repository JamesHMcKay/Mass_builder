#ifndef CALC_COUNTER_TERMS_H
#define CALC_COUNTER_TERMS_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"

using namespace utils;

class Calc_counter_terms
{
private:


  // definitions carried over from other file, remove soon
  std::map<std::string, Bases> full_basis;
  vector<string> full_basis_id;
  int np;
  
  Options options;
  string tag;
  string model;
  vector<string> masses_input,id_input;
  int nb;
  int nbr;
  vector<string> reduced_basis_id;
  std::map <std::string, Bases > reduced_basis;
  std::map <std::string, Bases > prod_basis;
  vector<std::string> prod_id;
  std::map <std::string, Bases_product > products_map;
    
  const char *ext = ".txt";
  string underscore = "_";
  string blank = "";

public:

  Calc_counter_terms(){}

  bool calc_counter_terms(Options options_in);
  std::string solve_1loop(std::string particle_name,vector<std::string> tags);
  
};

#endif