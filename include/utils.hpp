#ifndef UTILS_H
#define UTILS_H


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <ctime>

#include <unistd.h>

#include "bases.hpp" // to go in printer header later
#include "options.hpp"


using namespace std;

// need to divide utils up into actual utility like functions and printer functions


namespace utils
{

void user_input_guide();

void update_avail_diagrams(Options options);

const char * output_file_name(std::string model, std::string tag, std::string file);

void get_data(vector<std::string> &A,int &n,const char *filename);

void get_data(vector<std::string> &A, vector<std::string> &B,int &n,const char *filename, bool whole_line = false);

void get_data(vector<std::string> &A,vector<std::string> &B,vector<std::string> &C,int &n,const char *filename);

void print_math_header(ofstream &myfile);

void print_math_body(ofstream &file,Options options,string cwd);



bool check_done();

string part_name_simple(std::string particle_name_full);

vector<string> remove_duplicates(vector<string> input,string name);
vector<string> remove_duplicates(vector<string> input);

vector<char> remove_duplicates(vector<char> input,string name);
vector<char> remove_duplicates(vector<char> input);


string char_to_string(char c);

std::vector<std::string> extract_keys(std::map<std::string, Bases> const& input_map);

vector<int> find_string_lengths(vector<string> input);

void print_base(ofstream &myfile, Bases base, string id, string SEn);

void print_math_basis(std::map<std::string, Bases> base_map, ofstream &myfile , string target);

void print_math_products(std::map<std::string, Bases> base_map, ofstream &myfile, string target);


void ReplaceAll(std::string &input, const std::string& from, const std::string& to);

void print_doTSIL(ofstream &myfile,Bases base);

}

#endif