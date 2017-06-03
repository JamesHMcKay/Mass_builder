#ifndef GENERATE_CODE_H
#define GENERATE_CODE_H

// useful functions that are super specific to any area of the code


#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <cmath>
#include <ctime>

#include "options.hpp"
#include "templates.hpp"
#include "utils.hpp"
#include "cmake_variables.hpp"


using namespace std;
using namespace utils;


class  Generate_code
{
    time_t t = time(0);   // get time now to print into generated files
    struct tm * now = localtime( & t );
  
    string particle_name = "";
    string particle_name_reduced = "";
    string tag = "empty";
    
    const char* coeff_integrals_tmp = "/output/coeff_integrals_"; // vector containing file names
    const char* coeff_products_tmp = "/output/coeff_products_"; // vector containing file names
    const char* summation_tmp = "/output/summation_"; // vector containing file names
    
    const char *ext = ".txt";
    const char *cpp = ".cpp";
    string underscore = "_";
    string c_file_diagrams = "";
    string model;
  
    const char *file_diagrams;
    Options options;
  
    const char *file_masses;
    vector<string> masses_input,id_input;
    int na;
    vector<std::string> integrals;
  
    std::map<std::string, Bases> base_map;
    vector<std::string> couplings, relationships;
    int nc; // number of diagrams, length of diagram_number vector
    int nm; // number of diagrams, length of diagram_number vector
    int n_group; // number of diagrams in each source file
    vector<std::string> temp_vec; // TODO remove!
    vector<std::string> masses;
  
  
    // get list of available diagrams
    vector<std::string> tags;
    vector<std::string> particle_names,levels,particle_names_reduced;
    string level;
    int nd,np; // number of diagrams
    ofstream se_hpp;
    vector<std::string> particle_names_short_reduced;
    vector<int> subgrouplist;

    public:
    Generate_code(Options options) : options(options)
    {
      c_file_diagrams = options.input_list;
      file_diagrams = c_file_diagrams.c_str();
      model = options.model;
      get_data(particle_names,tags,levels, nd,file_diagrams);
      particle_names_reduced = particle_names;
      sort(particle_names_reduced.begin(),particle_names_reduced.end());
      particle_names_reduced.erase( unique( particle_names_reduced.begin(), particle_names_reduced.end() ), particle_names_reduced.end() );
      np = particle_names_reduced.size();
      cout << "np = " << np << endl;
      
      se_hpp.open ("include/self_energy.hpp");
      
      n_group = options.n_group;
      
    };
  
    void generate_code();
  
    void generate_particle_src(std::string particle,int);
  
    void decalare_var(ofstream &main_output);
    void decalare_var_tsil(ofstream &main_output);
  
    void generate_data_hpp();
  
    void generate_self_energy_hpp();
  
    pair<int,int> number_of_diagrams(std::string particle)
    {
      int num_diagrams = 0;
      for (int d = 0; d<nd;d++)
      {
        if (particle == particle_names[d])
        {
          num_diagrams = num_diagrams + 1;
          
          int diff = num_diagrams % n_group;
          if ((num_diagrams - diff)!=0)
          {
            subgrouplist[d] = (num_diagrams - diff)/n_group;
          }
          else
          {
            subgrouplist[d] = 0;
          }
          
        }
      }
      
      int modulo = num_diagrams % n_group;
      int groups = (num_diagrams - modulo)/n_group;
      
      pair<int,int> result;
      
      result.first = groups;
      result.second = modulo;
      
      return result;
    }
  
  
  
  
  
};
#endif