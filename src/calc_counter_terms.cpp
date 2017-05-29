/*
 Mass Builder 
 
 James McKay
 Feb 2017
 
 --- calc_counter_terms.cpp ---
 
 */

#include "calc_counter_terms.hpp"

#define RUN_ALL
//#define DEBUG
using namespace std;
using namespace utils;
string s_cwd_(getcwd(NULL,0));


std::string Calc_counter_terms::solve_1loop(std::string particle,vector<std::string> diagram)
{


  int nd = diagram.size();

  string dimension = "D";
  
  
  
  ofstream math_ct;
  math_ct.open ("output/math_ct.m");
  
  templates::print_math_header(math_ct);
  math_ct << "  SEtotal = 0 ;" << endl;
  
  for (int i=0;i<nd;i++)
  {
  
  math_ct<<"Get[\"" << s_cwd_ <<"/models/"<< options.model << "/output/math_data_" << particle << "_" << diagram[i] << "_1" << ".mx\"]\n";
 
  math_ct<<"SEtotal = SEnFinite + SEtotal;"<<endl;
  }
  //math_ct << "Print[\" -------------Total SE is ------------- \"]" <<endl;
  //math_ct << "Print[SEtotal]"<<endl;
  
 
  math_ct<<"SE = Coefficient[SEtotal,epsilon,-1]; \n";
  math_ct<<"SE = Simplify[SE /. epsilon->0];\n"; // some integrals come through as D = 4-epsilon so fix these
  
  
  // check for higher orders in 1/epsilon
  math_ct<<"SEhot = Coefficient[SEtotal,epsilon,-2] + Coefficient[SEtotal,epsilon,-3]; \n"
  << "Print[\" ------------- higher order terms in 1/epsilon ------------- \"]\n"
  << "Print[SEhot]\n"
  << "ct = FullSimplify[-SE*Pi^2/.Pair[Momentum[p], Momentum[p]]->p^2, {v == (2 mw/g2), g1 == (g2 STW/CTW), CTW^2 + STW^2 == 1,CTW == mw/mz}];\n"
  << "Print[\" ------------- 1-loop counter-term coupling is ------------- \"]\n"
  << "Print[ct]\n";
  math_ct.close();
 
#ifdef RUN_ALL
  system("chmod +x output/math_ct.m ");
  if(options.verbose) system("./output/math_ct.m");
  else system("./output/math_ct.m  >/dev/null ");
#endif

  
  std::string null;
  return null;
}





bool Calc_counter_terms::calc_counter_terms(Options options_in)
{
  
  bool success=0;
  options = options_in;
  
  // need to read in the list of available diagrams and then select the 1-loop ones for
  // adding up here
  
  // read in available diagrams

  const char *ext = ".txt";
  const char* file_diagrams_tmp = "models/";
  string c_file_diagrams = file_diagrams_tmp + options.model + "/output/avail_diagrams_" + ext;
  const char *file_diagrams = c_file_diagrams.c_str();

  vector<std::string> tags;
  vector<std::string> particle_names,levels;
  string level;
  int nd; // number of diagrams

  cout << "input list = " << file_diagrams << endl;

  get_data(particle_names,tags,levels, nd,file_diagrams);

  // one-loop corrections
  vector<std::string> tags_1;
  vector<std::string> particle_names_1,levels_1;
  
  for (int i=0;i<nd;i++)
  {
    if (levels[i]=="1" && particle_names[i]==options.particle)
    {
      particle_names_1.push_back(particle_names[i]);
      tags_1.push_back(tags[i]);
      levels_1.push_back(levels[i]);
    }
  
  }
  // send new array to Mathematica script which adds up
  // the amplitudes via another loop
  
  
  string particle_simple = part_name_simple(options.particle_1,options.particle_2);
  
  solve_1loop(particle_simple,tags_1);
  
  
  
  
  
     return success;
}
