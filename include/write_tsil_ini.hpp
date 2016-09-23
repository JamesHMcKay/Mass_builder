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

// this encodes a TSIL_EVALUATE call, so must have a definition of each mass

std::vector<Bases> integrals;
std::vector<string> masses;

public:

string x = "";
string y = "";
string z = "";
string u = "";
string v = "";

vector<bool> check_vec;
vector<int> location;

eval_obj() {}
  
eval_obj(string x, string y, string z, string u, string v) : x(x), y(y), z(z), u(u), v(v) {}

// create function which outputs all possible integrals this object is capable of evaluating

// since one of the above strings can be left blank, in that case the above function needs a list of possible masses


std::vector<Bases> get_integrals(std::vector<string> masses_input);
void add_integral(string type, string x, string y = "", string z = "", string u = "", string v = "");

std::vector<bool> get_check_vec(vector<string> names, std::map<std::string, Bases> base_map,std::vector<string> masses_input);


};






class Print_dotsil
{

  vector<string> names;
  std::map<std::string, Bases> base_map;
  vector<string> masses;
  vector<vector<bool>> V_check_vec;
  std::vector<eval_obj> eval_vec;

public:
  
  Print_dotsil() {}
  Print_dotsil(vector<string> names, std::map<std::string, Bases> base_map) : names(names) , base_map(base_map) {}
  
  void print_to_file(ofstream &myfile);
  void sort_integrals();
  
  void get_poss_eval(Bases base);
  void add_eval_obj( string x, string y, string z, string u, string v);
  void make_sets();
  
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