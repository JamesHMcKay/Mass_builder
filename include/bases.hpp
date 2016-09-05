#ifndef BASES_H
#define BASES_H

#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>


using namespace std;

class Bases
{

public:

string type = "";

string e1 = " ", e2 = " ", e3 = " ", e4 = " ", e5 = " ";

string coefficient = "";

string short_name = "";

Bases() {}

Bases(string type, string e1, string e2, string e3, string e4, string e5) : type(type) , e1(e1), e2(e2), e3(e3), e4(e4), e5(e5) {}

Bases(string type, string e1, string e2, string e3, string e4) : type(type) , e1(e1), e2(e2), e3(e3),e4(e4) {}

Bases(string type, string e1, string e2, string e3) : type(type) , e1(e1), e2(e2), e3(e3) {}

Bases(string type, string e1, string e2) : type(type) , e1(e1), e2(e2) {}

Bases(string type, string e1) : type(type) , e1(e1) {}


void set_coeff(string _coefficient){coefficient = _coefficient;}

string get_coeff(){return coefficient;}

};

vector<int> find_string_lengths(vector<string> input);

string get_id(std::vector<string> &masses, std::vector<string> &identifiers, string mass);

void set_id(std::vector<string> &masses_input, std::vector<string> &identifiers_input);

string get_short_name(Bases basis, std::vector<string> &masses, std::vector<string> &identifiers);

std::map<std::string, Bases> set_bases(std::vector<string> masses, std::vector<string> &identifiers);

std::map <std::string, Bases > remove_zeros(std::map <std::string, Bases > base_map, std::vector<std::string> base_names);

void format_coeff(std::map <std::string, Bases > &base_map, std::vector<std::string> bases_names,std::vector<std::string> masses,std::vector<std::string> id);

std::map <std::string, Bases > products_container(vector<string> bases_names);



#endif