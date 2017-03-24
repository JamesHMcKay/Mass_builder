/*
 Mass Builder - the missing link in automated two-loop self energy calculations
 Please refer to the documentation for details or the readme.txt for simple run instructions
 
 James McKay
 
 FeynArts (hep-ph/0012260) is used to generate the two-loop amplitudes
 FeynCalc (ArXiv:1601.01167) is used to reduce the amplitudes
 TARCER (hep-ph/9801383) is used to reduce the resulting amplitudes to basis integrals
 TSIL (hep-ph/0501132) is used to evaluate the basis integrals
 
 */
#include "data.hpp"
#include "calc_amplitudes.hpp"
#include "calc_counter_terms.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"
//#include "supplements.hpp"
//#include "write_tsil_ini.hpp"

using namespace std;
//using namespace supplementary_code;
using namespace utils;



void run_mass_builder_mode_1a(Options options)
{
  Calc_amplitudes ca;
  std::string particles [1000];
  std::string types [1000];
  std::string diagrams [1000]; int i=0;
  std::string model = options.model;
  const char *file_diagrams;
  if (options.input_list==""){
    const char *ext = ".txt";
    const char* file_diagrams_tmp = "models/";
    string c_file_diagrams = file_diagrams_tmp + model + "/diagrams" + ext;
    file_diagrams = c_file_diagrams.c_str();
  }
  else
  {
    file_diagrams = options.input_list.c_str();
  }
  std::ifstream input(file_diagrams);
  std::string line;
  while(getline(input, line))
  {
    if (!line.length() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    iss>> particles[i] >> diagrams[i] >> types[i];
    i=i+1;
  }
  for (int k=0;k<i;k++)
  {
    options.particle = particles[k];
    options.diagram = diagrams[k];
    options.set_type(types[k]);
    options.model = model;
    ca.calc_diagram(options);
  }
}

void run_mass_builder_mode_1b(Options options)
{
  Calc_amplitudes ca;
  ca.calc_diagram(options);
}

void run_mass_builder_mode_2(Options options)
{
  Calc_counter_terms ct;
  ct.calc_counter_terms(options);
}

void run_mass_builder_mode_6(Options options)
{
  Data data(options);
  Self_energy se;
  
  se.run_tsil(data);
  for (unsigned int i = 0;i < data.avail_part.size();i++)
  {
    cout << "One loop self energy of particle " << data.avail_part[i] << " = " << data.SE_1[data.avail_part[i]] << endl;
    cout << "Two loop self energy of particle " << data.avail_part[i] << " = " << data.SE_2[data.avail_part[i]] << endl;
  }
  
  
}


void run_mass_builder_mode_7(Options options)
{
  
  // iterate to determine physical mass of spin 0 particle
  Data data(options);
  for (unsigned int i = 0;i < data.avail_part.size();i++)
  {
    
    
    //double M_tree = pow(data.M_tree[data.avail_part[i]],2);
    double M_tree = pow(data.M_tree[data.avail_part[i]],2);
    
   
   // cout << "data.M_tree[data.avail_part[i]] = " << data.M_tree[data.avail_part[i]] << endl;
    
    double diff = 10;
    data.P = data.M_tree[data.avail_part[i]];
    double mass = 0;
    
    while (diff>0.0001)
    {
      Self_energy se;
      
      se.run_tsil(data);
      
      double self_energy = data.SE_1[data.avail_part[i]];
      double mass_sqr = M_tree - self_energy;
      
      mass = pow(abs(mass_sqr),0.5);
      
      diff = abs(mass - data.P);
     // cout << "diff = " << diff << endl;
     // cout << "self energy = " << self_energy << endl;
     // cout << "M_tree = " << pow(M_tree,0.5) << endl;
     // cout << "M_physical = " << mass << endl;
      data.P = mass;
    }
    
    cout << "physical one-loop mass of " << data.avail_part[i] << " = " << mass << endl;
    
  }
  
}



int main(int argc, char *argv[])
{
  
  // testing
  
  //Print_dotsil dotsil;
  //dotsil.sort_integrals();
  
  
  
  
  User_input user(argc,argv);
  
  user.user_interface();
  
  Options options = user.options;
  
  // read options and work through possibilities for each run mode and check requirements are meant
  
  if (options.model=="" && (options.run_mode < 6)){ cout << "no model specified" << endl; return 0;}
  
  if (options.run_mode == 1)
  {
    if ((options.particle == "") || (options.diagram == ""))
    {
      run_mass_builder_mode_1a(options);
    }
    else
    {
      if ((options.input_list==""))
      {
        if ((options.particle == "") || (options.diagram == "")) { cout << "no valid input selected" << endl;}
        else { run_mass_builder_mode_1b(options);}
      }
    }
  }
  
  if (options.run_mode == 2)
  {
    if ((options.particle == "")) { cout << "please enter a particle" << endl; return 0;}
    else { run_mass_builder_mode_2(options);}
  }
  
  if (options.run_mode == 4 )
  {
    if (options.input_list == "")
    {
      options.input_list = "models/"+ options.model + "/output/avail_diagrams_.txt";
      sort_avail_diagrams(options);
    }
    if (options.model == "") { cout << "please specify a model" << endl; return 0;}
    Generate_code gen_code(options);
    gen_code.generate_code();
  }
  
  if (options.run_mode == 5 )
  {
    if (options.model == "" || options.particle == "") { cout << "please specify a model and particle, at least one is missing" << endl; return 0;}
    Calc_amplitudes ca;
    ca.generate_figures(options);
  }
  
  if (options.run_mode == 6 )
  {
    if (options.input_list == "" ) { cout << "missing input data" << endl; return 0;}
    run_mass_builder_mode_6(options);
  }
  
  if (options.run_mode == 7 )
  {
    if (options.input_list == "" ) { cout << "missing input data" << endl; return 0;}
    run_mass_builder_mode_7(options);
  }
  
}
