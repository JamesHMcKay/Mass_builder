#ifndef COMPUTE_AMP_H
#define COMPUTE_AMP_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"
#include "templates.hpp"
#include "cmake_variables.hpp"

#include "wstp.h"
#include<math.h>
#include <stdio.h>

using namespace utils;

class Compute_amp
{
private:
  
  std::map<std::string, Bases> full_basis;
  vector<string> full_basis_id;
  int np;
  bool multi_particle;
  bool seperate_amp = false;
  ofstream debug_out;
  ofstream log_out;
  
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
  
  void *pHandle;
  WSLINK link;
  
public:
  
  Compute_amp(Options options) : options(options) {}
  
  bool calc_diagram();
  
  
  void create_wstp_link()
  {
    
    int WSerrno;
    WSENV WSenv = WSInitialize(0);
    if(WSenv == (WSENV)0)
    {
      cout << "Unable to initialize WSTP environment" << endl;
      
      
    }
    else
    {
      cout << "The environment is initialized successfully..." << endl;
    }
    
    
    std::stringstream WSTPflags;
    
    // This opens a WSTP connection
    #ifdef __APPLE__
      WSTPflags << "-linkname " << MATHEMATICA_KERNEL << " -mathlink";
    #else
      WSTPflags << "-linkname math -mathlink";
    #endif
    
    pHandle = WSOpenString(WSenv, WSTPflags.str().c_str(), &WSerrno);
    
    link = (WSLINK)pHandle;
    
    if(link == (WSLINK)0 || WSerrno != WSEOK)
    {
      cout << "Unable to create link to the Kernel" << endl;
      WSNewPacket(link);
    }
    else
    {
      cout << "WSTP link started" << endl;
    }
    
    
  }
  
  // Load required Mathematica libraries
  void load_libraries()
  {
		string cwd = getcwd(NULL,0);
    
    std::string input ="";
    
    input = "$LoadAddOns = {\"FeynHelpers\"};";
    input+= "$LoadTARCER = True;";
    input+= "$LoadFeynArts = True;";
    input+= "<< FeynCalc/FeynCalc.m;";
    input+= "AppendTo[$Path, \"" + cwd + "/src/\"];";
    input+= "<< MassBuilder.m;";
		
    if (options.verbose)
    {
			ReplaceAll(input,";","; \n ");
			debug_out << input << endl;
			debug_out << "   " << endl;
			debug_out << "(* ::Section:: *)" << endl;
		}
		else
		{
	    WSNewPacket(link);
	    WSPutFunction(link, "ToExpression", 1);
	    
	    WSPutString(link, input.c_str());
	    
	    wait_for_packet();
	   }
	   if (!options.verbose)
	   {
		   // send to log
		   ReplaceAll(input,";","; \n ");
			 log_out << input << endl;
			 log_out << "   " << endl;
			 log_out << "(* ::Section:: *)" << endl;
		}
		
  }
  
  // Wait to receive a packet from the kernel
  void wait_for_packet()
  {
    int pkt;
    while( (pkt = WSNextPacket(link), pkt) && pkt != RETURNPKT)
    {
      WSNewPacket(link);
      if (WSError(link))
      {
        cout << "Error reading packet from WSTP" << endl;
      }
    }
  }
  
  // Send a string to be evaluated in Mathematica via WSTP
  void send_to_math(std::string &input)
  {
    if (options.verbose)
    {
			ReplaceAll(input,";","; \n ");
      debug_out << input << endl;
    }
    else
	  {
	    WSNewPacket(link);
	    WSPutFunction(link, "ToExpression", 1);
	    WSPutString(link, input.c_str());
	    
	    wait_for_packet();
		}
		if (!options.verbose)
		{
			ReplaceAll(input,";","; \n ");
      log_out << input << endl;
		}
    input = "";
  }
  
  void open_log_files()
  {
	  if (options.verbose)
	  {
			debug_out.open ("debug.m");
			debug_out << "(* ::Package:: *)" << endl;
			debug_out << "Quit[]" << endl;
			debug_out << "(* ::Section:: *)" << endl;
		}
		else
		{
			log_out.open ("output/log.m");
			log_out << "(* ::Package:: *)" << endl;
			log_out << "Quit[]" << endl;
			log_out << "(* ::Section:: *)" << endl;
		}
	}
  
  
  
  // Generate all loop diagrams for specified particle
  void generate_figures();
  
  // Compute tree-level counter-term coupling (print result to terminal)
  void calc_counter_terms();
  
  // solve the 1-loop order linear equation to determine counter-term coupling
  void solve_1loop(std::string particle,vector<std::string> diagram);
  
  // temporary functions for making a quick list of 2 loop diagrams
  void make_2loop_list(std::string particle,vector<std::string> diagram,vector<std::string> level);
    
  // Main function to print vertices and Feynman rules to LaTeX ready file
  void print_vertices();
  
  // Generate tex file
  void print_vertex_tex(ofstream &myfile, std::string particle_1,std::string particle_2,std::string particle_3,std::string particle_4,Options options,int number);
  
  // print all relevant four vertices to file
  void print_4_vertex(std::string &input, std::string particle_1,std::string particle_2,std::string particle_3,std::string particle_4,Options options);
  
  // print all relevant three vertices to file
  void print_3_vertex(std::string &input, std::string particle_1,std::string particle_2,std::string particle_3,Options options);

};

#endif
