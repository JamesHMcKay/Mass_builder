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

class eval_obj
{
  // this object encodes a TSIL_EVALUATE call

  std::vector<Bases> integrals;
  vector<bool> check_vec; // what required basis integrals does this statement fufill
  bool check_vec_evaluated = false;

public:

  string x = "";
  string y = "";
  string z = "";
  string u = "";
  string v = "";

  
  
  
  vector<int> location;   // what is the location of this integral within this object (will be required to hold the string passed to TSIL_EVALUATE)
  std::vector<string> eval_string;



  eval_obj() {}
    
  eval_obj(string x, string y, string z, string u, string v) : x(x), y(y), z(z), u(u), v(v) {}

  std::vector<Bases> get_integrals();
  void add_integral(string type, string short_name, string x, string y = " ", string z = " ", string u = " ", string v = " ");

  std::vector<bool> get_check_vec(vector<string> names, std::map<std::string, Bases> base_map);
  void swap_tsil_to_tarcer_V(string &a,string &b,string &c, string &d);
};



class Print_dotsil
{

  vector<string> names;
  std::map<std::string, Bases> base_map;
  vector<string> masses;
  vector<vector<bool>> V_check_vec;
  std::vector<eval_obj> eval_vec;
  int eval_count = 0;

public:
  
  Print_dotsil() {}
  Print_dotsil(vector<string> names, std::map<std::string, Bases> base_map) : names(names) , base_map(base_map) {}
  
  void print_to_file(ofstream &myfile);
  void sort_integrals();
  
  void get_poss_eval(Bases base);
  void add_eval_obj( string x, string y, string z, string u, string v);
  void make_sets();
  void print_total(std::vector<int>);
  string coeff(string type);
  
  void print_eval_obj(ofstream &myfile,eval_obj &eo, vector<int> &total);
  
  vector<int> get_duplicates();
  
  
  // indepdently determine a unique list of all masses involved in this problem
  void get_masses()
  {
    for (unsigned int i=0;i<names.size();i++)
    {
      if (base_map[names[i]].e1!=" "){masses.push_back(base_map[names[i]].e1);}
      if (base_map[names[i]].e2!=" "){masses.push_back(base_map[names[i]].e2);}
      if (base_map[names[i]].e3!=" "){masses.push_back(base_map[names[i]].e3);}
      if (base_map[names[i]].e4!=" "){masses.push_back(base_map[names[i]].e4);}
      if (base_map[names[i]].e5!=" "){masses.push_back(base_map[names[i]].e5);}
    }
    
    sort(masses.begin(),masses.end());
    masses.erase( unique( masses.begin(), masses.end() ), masses.end() );
  };
  
  
};

#endif