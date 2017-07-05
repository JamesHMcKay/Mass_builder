#ifndef UTILS_H
#define UTILS_H

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

#include "bases.hpp"
#include "options.hpp"


using namespace std;

namespace utils
{
  
  int get_loop_order(string type);
  
  string get_cwd();
  
  const char * add_mpi_ext(std::string name, int process, std::string ext);
  
  string char_to_string(char c);
  
  vector<string> extract_keys(map<string, Bases> const& input_map);
  
  vector<int> find_string_lengths(vector<string> input);
  
  bool check_done(std::string &result, int mpi_process);
  
  bool check_done_quiet(int mpi_process);
  
  string part_name_simple(string particle_name_full);
  
  string part_name_simple(std::string particle_name_full_1,std::string particle_name_full_2);
  
  vector<string> remove_duplicates(vector<string> input,string name);
  vector<string> remove_duplicates(vector<string> input);
  
  vector<char> remove_duplicates(vector<char> input,string name);
  vector<char> remove_duplicates(vector<char> input);
  
  void update_avail_diagrams(Options options);
  
  string make_tag(Options &options);
  
  bool check_if_available(Options options);
  
  void sort_avail_diagrams(Options options);
  
  void ReplaceAll(string &input, const string& from, const string& to);
  
  const char * output_file_name(string model, string tag, string file);
  
  void get_data(vector<string> &A,int &n,const char *filename);
  
  void get_data(vector<string> &A, vector<string> &B,int &n,const char *filename, bool whole_line = false);
  
  void get_data(vector<string> &A,vector<string> &B,vector<string> &C,int &n,string filename);
  
  void get_data(vector<std::string> &A,vector<std::string> &B,vector<std::string> &C,vector<std::string> &D,int &n, string filename);
  
  void assign_FCGV(ofstream &file,Options options);
  
  void assign_FCGV(std::string &file,Options options);
  
  void assign_variables(ofstream &file,Options options);
  
  void assign_variables(std::string &file,Options options);
  
  void get_saved_amplitude(std::string &input, Options options);
  
  void print_math_body(std::string &input,Options options,string cwd,std::vector<std::string> masses);
  
  void print_tarcer_recurse(std::string &input);
  
  void print_base(ofstream &myfile, Bases base, string id, string SEn,string D);
  
  void get_finite_amplitude(std::string &input, Options options);
  
  void print_math_basis(map<string, Bases> base_map, ofstream &myfile , string target, string D);
  
  void print_math_products(map<string, Bases> base_map, ofstream &myfile, string target, string D);
  
  void print_doTSIL(ofstream &myfile,Bases base);
  //void print_doTSIL_cout(Bases base);
  
  void print_finite_base(ofstream &myfile, Bases base, string id);
  
  void print_finite_basis(std::map<std::string, Bases> base_map, ofstream &myfile);
  
  void timestamp();
  
  void print_diagram_info(Options options);
  
}

#endif