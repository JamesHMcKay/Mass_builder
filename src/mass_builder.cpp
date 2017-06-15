/*
 Mass Builder
 
 James McKay
 Aug 2016 - June 2017
 
 --- mass_builder.cpp ---
 
 This file manages the execution of all possible tasks within the Mass Builder program
 by passing the user input to the user interface routine and then interpreting the resultant
 Options class that is returned.
 
 In this file we manage the MPI execution of a many diagram calcuation.
 
*/

# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <mpi.h>

#include "data.hpp"
#include "generate_code.hpp"
#include "self_energy.hpp"
#include "compute_amp.hpp"

using namespace std;
using namespace utils;

// subroutine used by run_mass_builder_mode_1a to manage MPI tasks
bool check_task_done(int number_of_processes, int task_number)
{
  for (int i = 0; i<number_of_processes ; i++)
  {
    
    const char* file_tmp = "output/tasks_";
    string c_file = file_tmp + std::to_string(i) + ".txt";
    const char *file = c_file.c_str();
    vector<string> completed_tasks;
    int length;
    get_data(completed_tasks, length, file );
    for (int j = 0; j<length ; j++)
    {
      if (completed_tasks[j] == std::to_string(task_number))
      {
        return true;
      }
    }
  }
  return false;
}

// subroutine used by run_mass_builder_mode_1a to manage MPI tasks
void inform_task_started(int processes_number, int task)
{
  const char* file_tmp = "output/tasks_";
  string c_file = file_tmp + std::to_string(processes_number) + ".txt";
  const char *file = c_file.c_str();
  vector<string> completed_tasks;
  int length;
  get_data(completed_tasks, length, file );
  
  ofstream outfile;
  outfile.open (file);
  
  for (int j = 0; j<length ; j++)
  {
    outfile << completed_tasks[j] << endl;
  }
  
  outfile << task;
  
  outfile.close();
}

// compute a list of diagrams
void run_mass_builder_mode_1a(Options options,int argc, char *argv[])
{
  int id;
  int p; // number of processes
  double wtime = 0;
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
  
  
  //  Initialize MPI
  MPI_Init ( &argc, &argv );
  //  Get the number of processes.
  MPI_Comm_size ( MPI_COMM_WORLD, &p );
  //  Get the individual process ID
  MPI_Comm_rank ( MPI_COMM_WORLD, &id );
  
  //  Process 0 prints initial message
  if ( id == 0 )
  {
    timestamp ( );
    cout << "\n";
    cout << "Mass Builder is running\n";
    cout << "\n";
    cout << "  The number of processes is " << p << "\n";
    cout << "\n";
  }
  
  if ( id == 0 )
  {
    wtime = MPI_Wtime ( );
  }
  
  bool done = false;
  int task = id;
  int original_task = id;
  
  if ( p > i )
  {
    cout << "There are more processes running than tasks to complete, please rerun with less processes or request more diagrams." << endl;
  }
  else
  {
    while ( done == false )
    {
      // if task not yet started or completed by another process then go ahead
      if ( !check_task_done(p,task) )
      {
        // write to file to inform other processes that this task has started
        inform_task_started(id, task);
        
        // execute calculation
        Compute_amp ca(options);
        options.particle_1 = particles[task];
        options.particle_2 = particles[task];
        options.particle = particles[task];
        options.diagram = diagrams[task];
        options.mpi_process = id;
        options.set_type(types[task]);
        ca.calc_diagram();
      }
      
      task = task + p;
      if ( !( task < i ) )
      {
        // reset num to be offset from original task id and start again
        task = (task + 1) % (p);
        
        if (task == original_task)
        {
          done = true;
          cout << "process " << id << " has completed all tasks" << endl;
        }
      }
    }
  }
  if ( id == 0 )
  {
    wtime = MPI_Wtime ( ) - wtime;
    cout << "  Elapsed wall clock time = " << wtime << " seconds.\n";
  }
  
  //  Terminate MPI
  MPI_Finalize ( );
  
  if ( id == 0 )
  {
    cout << "\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
    system("rm output/tasks*.txt >/dev/null");
  }
  
  
}

// compute one specified diagram
void run_mass_builder_mode_1b(Options options)
{
  //Calc_amplitudes ca; // alternative method
  Compute_amp ca(options);
  ca.calc_diagram();
}

// Compute tree-level counter-term coupling
void run_mass_builder_mode_2(Options options)
{
  Compute_amp ct(options);
  ct.calc_counter_terms();
}


// Evaluate self energies of all available particles
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


// iterate to determine physical mass of spin 0 particle
void run_mass_builder_mode_7(Options options)
{
  Data data(options);
  for (unsigned int i = 0;i < data.avail_part.size();i++)
  {
    double M_tree = pow(data.M_tree[data.avail_part[i]],2);
    
    double diff = 10;
    data.P = data.M_tree[data.avail_part[i]];
    double mass = 0;
    while (diff>0.0001)
    {
      Self_energy se;
      se.run_tsil(data);
      
      double self_energy = data.SE_1[data.avail_part[i]];
      double mass_sqr = M_tree + self_energy;
      
      mass = pow(abs(mass_sqr),0.5);
      
      diff = abs(mass - data.P);
      
      data.P = mass;
    }
    cout << "physical one-loop mass of " << data.avail_part[i] << " = " << mass << endl;
  }
}

int main(int argc, char *argv[])
{
  // pass input to user interface routine
  User_input user(argc,argv);
  
  user.user_interface();
  
  // set options class with user input
  Options options = user.options;
  
  // read options and work through possibilities for each run mode and check requirements are met
  
  if (options.model=="" && (options.run_mode < 6)){ cout << "no model specified" << endl; return 0;}
  
  // compute amplitudes using Mathematica tools
  if (options.run_mode == 1)
  {
    if ((options.particle == "") || (options.diagram == ""))
    {
      run_mass_builder_mode_1a(options,argc, argv);
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
  
  // compute tree-level counter-term coupling
  if (options.run_mode == 2)
  {
    if ((options.particle == "")) { cout << "please enter a particle" << endl; return 0;}
    else { run_mass_builder_mode_2(options);}
  }
  
  // generate C++ interface to TSIL
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
  
  // create pdf of FeynArts generated Feynman diagrams
  if (options.run_mode == 5 )
  {
    if (options.model == "" || options.particle == "") { cout << "please specify a model and particle, at least one is missing" << endl; return 0;}
    Compute_amp ca(options);
    ca.generate_figures();
  }
  
  // Evaluate self energies via TSIL interface
  if (options.run_mode == 6 )
  {
    if (options.input_list == "" ) { cout << "missing input data" << endl; return 0;}
    run_mass_builder_mode_6(options);
  }
  
  // Evaluate pole mass for a spin 0 field by iteration
  if (options.run_mode == 7 )
  {
    if (options.input_list == "" ) { cout << "missing input data" << endl; return 0;}
    run_mass_builder_mode_7(options);
  }
  
  // Print Feynman rules for specified vertices to LaTeX format
  if (options.run_mode == 8 )
  {
    if (options.model == "" ) { cout << "please specify a model to work with" << endl; return 0;}
    Compute_amp ca(options);
    ca.print_vertices();
  }
  
  return 0;
}
