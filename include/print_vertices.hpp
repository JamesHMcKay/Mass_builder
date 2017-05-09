#ifndef PRINT_VERTICES_H
#define PRINT_VERTICES_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"

using namespace utils;


void print_vertices(Options options_in);

void print_vertex(ofstream &myfile, std::string particle_1,std::string particle_2,std::string particle_3,Options options);


#endif