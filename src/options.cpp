/*
 Mass Builder
 
 James McKay
 Sep 2016
 
--- options.cpp ---
 
 The functions defined in this file deal read the user input from
 the command line and set up the Options class which is passed
 throughout the program
*/

#include "options.hpp"
#include "utils.hpp"
#include "templates.hpp"

void Options::set_type(std::string type)
{
  counter_terms = false;
  string loop_order_str = utils::char_to_string(type[0]);
  stringstream convert(loop_order_str);
  convert >> loop_order;
  string c_str = utils::char_to_string(type[1]);
  if (c_str == "c")
  {
    counter_terms = true;
  }
}

void User_input::call_user_guide()
{
  templates::user_input_guide();
}

void User_input::user_interface()
{
  string io;
  if (find_string("-c")){ options.counter_terms = true;}
  if (find_string("-v")){ options.verbose = true;}
  if (find_string("-o")){ options.optimise = true;}
  if (find_string("-w")){ options.detailed_output = true;}
  if (find_string("-a")){ options.run_mode = 1;}
  if (find_string("-b")){ options.run_mode = 2;}
  if (find_string("-s")){ options.n_group = 2;}
  if (find_string("-l"))
  {
    string input = "loop order";
    if (find_and_read_string("-l",input))
    {
      if (input == "2")
      {
				options.loop_order = 2;
			}
      else if (input == "1")
      {
				options.loop_order = 1; 
				options.onshell = false;
			}
			// not actually supported -- used as a temporary option
			else if (input == "3")
      {
				options.loop_order = 3; 
				options.onshell = false;
			}			
      else 
      {
				cout <<"This loop order is not supported please enter 1 or 2"<<endl;
			}
    }
  }

  if (find_string("-m"))
  {
    string input = "model name";
    if (find_and_read_string("-m",input)){options.model = input;}
    
    if (options.model.find("MSSM") != std::string::npos)
    {
      options.restrictions = "WinoLimit,WinoCouplings";
    }
    if (options.model == "VDM" || options.model == "QED"|| options.model == "QED2")
    {
      options.use_lorentz = false;
    }
    if (options.model == "MDM_triplet")
    {
      options.restrictions = "noFermions";
    }    

  }

  if (find_string("-p"))
  {
    string input = "particle name";
    if (find_and_read_string("-p",input))
    {
      options.particle = input;
      options.particle_1 = input;
      options.particle_2 = input;
    }
  }
  
  if (find_string("-q"))
  {
    string input = "secondary particle name";
    if (find_and_read_string("-q",input))
    {
      options.particle_2 = input;
    }
  }

  if (find_string("-d"))
  {
    string input = "diagram number";
    if (find_and_read_string("-d",input)){options.diagram = input;}
  }

  if (find_string("-i"))
  {
    string input = "a diagram list";
    if (find_and_read_string("-i",input)){options.input_list = input;}
  }

  if (find_string("-r"))
  {
    string input = "restrictions";
    if (find_and_read_string("-r",input)){options.restrictions = input;}
  }


  if (find_string("-g"))
  {
    string input = "generate code";
    if (find_and_read_string("-g",input))
    {
      options.run_mode = 4;
    }
  }
  
  
  if (find_string("-n"))
  {
    string input = "number of final states";
    if (find_and_read_string("-n",input))
    {
      if (input == "2"){options.n_final_states = 2;}
    }
  }
  

  if (find_string("-k"))
  {
    string input = "epsilon order";
    if (find_and_read_string("-k",input)){options.epsilon_order = input;}
  }

  if (find_string("-f")){ options.run_mode = 5;}
  if (find_string("-e")){ options.run_mode = 6;}
  if (find_string("-0")){ options.run_mode = 7;}
  if (find_string("-x")){ options.run_mode = 8;}
  if (find_string("-y")){ options.run_mode = 9;}
  
  if (find_string("-z")){ options.use_lorentz = false;}
  
  if (find_string("reuse")){ options.force = false;}
  
  if (find_string("onshell")){ options.onshell = false;}
  
  if (find_string("expole")){ options.exclude_photon_pole = false;}
  
  if (find_string("gauge")){ options.general_gauge = true; options.fire = false;}
  
  if (find_string("-t")){ options.fire = false;}
  
  if (options.particle_1!=options.particle_2)
  {
    options.particle = options.particle_1 + "_" + options.particle_2;
  }
  

  #ifdef DEBUG
  options.print_options();
  #endif
  
}


