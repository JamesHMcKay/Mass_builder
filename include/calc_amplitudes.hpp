#ifndef CALC_AMPLITUDES_H
#define CALC_AMPLITUDES_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"
#include "templates.hpp"

using namespace utils;

class Calc_amplitudes
{
private:

  std::map<std::string, Bases> full_basis;
  vector<string> full_basis_id;
  int np;
  bool multi_particle;
  
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

  Calc_amplitudes(){}

  void compute_amp(string prevb,string dimension);
  void make_finite_amp(bool counter_terms);
  void make_full_trial(string prevb,string dimension,bool cform);
  bool calc_diagram(Options options_in);
  void generate_figures(Options options_in);
  void initial_trial(string dimension);
  void second_initial_trial(string prevb,string dimension);
};

#endif