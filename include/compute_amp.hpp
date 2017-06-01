#ifndef COMPUTE_AMP_H
#define COMPUTE_AMP_H

#include "utils.hpp"
#include "bases.hpp"
#include "options.hpp"
#include "templates.hpp"

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
  
  Compute_amp(){}
  
  bool calc_diagram(Options options_in);
  
  
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
    
    
    /*MATH_KERNEL_PATH */  WSTPflags << "-linkname " << "/Applications/Mathematica.app/Contents/MacOS/MathKernel" << " -mathlink";
    
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
  
  void load_libraries()
  {
    
    WSNewPacket(link);
    WSPutFunction(link, "ToExpression", 1);
    
    string cwd = getcwd(NULL,0);
    
    std::string input;
    
    input = "$LoadTARCER = True;";
    input+= "$LoadFeynArts = True;";
    input+= "<< FeynCalc/FeynCalc.m;";
    input+= "AppendTo[$Path, \"" + cwd + "/src/\"];";
    input+= "<< MassBuilder.m;";
    
    WSPutString(link, input.c_str());
    
    wait_for_packet();
  }
  
  void wait_for_packet()
  {
  /* Wait to receive a packet from the kernel */
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
  
  void send_to_math(std::string &input)
  {
    WSNewPacket(link);
    WSPutFunction(link, "ToExpression", 1);
    WSPutString(link, input.c_str());
    wait_for_packet();
    input = "";
  }
  
  
  void generate_figures(Options options_in);
  
  void calc_counter_terms(Options options_in);
  
  void solve_1loop(std::string particle,vector<std::string> diagram);
  
  
  
};

#endif
