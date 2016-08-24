

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


#define RUN_ALL

using namespace std;


void draw_all_diagrams(std::string particle, string model)
{

  bool verbose=1;




  
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );
  
  
  
  
  string s_cwd(getcwd(NULL,0));
  
  
  
  ofstream myfile;
  myfile.open ("make_figures.m");


  myfile<< "#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script\n"
  << "(* ---------------------------------------------------- *)\n"
  << "(* This file has been automatically generated by stage_1.cpp, on "<< now->tm_mday << '-'
  << (now->tm_mon + 1) << '-'<< (now->tm_year + 1900) <<", do not edit *)\n"
  << "(* ---------------------------------------------------- *)\n"
  << "$LoadPhi = True;\n"
  << "$LoadTARCER = True;\n"
  <<"$LoadFeynArts = True;\n"
  << "<< FeynCalc/FeynCalc.m\n"
  <<"dm[mu_] := DiracMatrix[mu, Dimension -> D]\n"
  <<"dm[5] := DiracMatrix[5]\n"
  <<"ds[p_] := DiracSlash[p]\n"
  <<"SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True];\n"
  <<"SetOptions[DiracTrace, DiracTraceEvaluate -> True];\n"
  <<"$GenericMixing = True;\n"
  <<"t12 = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> Internal];\n"
  <<"alldiags = InsertFields[t12, {"<<particle<<"} -> {"<<particle<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<s_cwd<<"/models/"<<model<<"/"<<model<<"\"];\n"
  <<"Export[\""<<s_cwd<<"/FA_diagrams/all_diagrams_"<<particle<<".pdf\",Paint[alldiags]];\n"  // print the FA diagram to pdf in local directory
  <<endl;

  #ifdef RUN_ALL
  system("chmod +x make_figures.m ");
  if (verbose) system("./make_figures.m");
  else system("./make_figures.m  >/dev/null");
  #endif
  
  
  
}




void draw_diagrams(vector<std::string> particles, vector<std::string> diagrams, int nd,string model)
{

  bool verbose=1;




  
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );
  
  
  
  
  string s_cwd(getcwd(NULL,0));
  
  
  
  
  
  ofstream myfile;
  myfile.open ("make_figures.m");


  myfile<< "#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script\n"
  << "(* ---------------------------------------------------- *)\n"
  << "(* This file has been automatically generated by stage_1.cpp, on "<< now->tm_mday << '-'
  << (now->tm_mon + 1) << '-'<< (now->tm_year + 1900) <<", do not edit *)\n"
  << "(* ---------------------------------------------------- *)\n"
  << "$LoadPhi = True;\n"
  << "$LoadTARCER = True;\n"
  <<"$LoadFeynArts = True;\n"
  << "<< FeynCalc/FeynCalc.m\n"
  <<"dm[mu_] := DiracMatrix[mu, Dimension -> D]\n"
  <<"dm[5] := DiracMatrix[5]\n"
  <<"ds[p_] := DiracSlash[p]\n"
  <<"SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True];\n"
  <<"SetOptions[DiracTrace, DiracTraceEvaluate -> True];\n"
  <<"$GenericMixing = True;\n"
  <<" (*MDM_tripletEWSB*)\n"
  <<"t12 = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> Internal];\n"
  //<<"chi0 = InsertFields[t12, {F[6]} -> {F[6]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  //<<"chi1 = InsertFields[t12, {F[5]} -> {F[5]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  <<endl;
  
  
  vector<std::string> particle_names_short = particles;


  sort(particle_names_short.begin(),particle_names_short.end());
  particle_names_short.erase( unique( particle_names_short.begin(), particle_names_short.end() ), particle_names_short.end() );
  bool first = 1;
  for (int i=0;i<particle_names_short.size();i++)
  {
  string particle_name_tmp = particle_names_short[i];
  myfile<<particle_name_tmp<<" = InsertFields[t12, {"<<particle_name_tmp<<"} -> {"<<particle_name_tmp<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<s_cwd<<"/models/"<<model<<"/"<<model<<"\"];\n";
  myfile<< "subdiags" << particle_name_tmp <<" =   DiagramExtract["<<particle_name_tmp;
  for (int d = 0; d<nd;d++)
  {
  if (particles[d] == particle_name_tmp)
  {
  //if (first) myfile <<", ( "<<diagrams[d];
  myfile <<", "<<diagrams[d];
  first =0;
  }
  }
  myfile <<"]\n";
  first = 1;
  myfile <<"Export[\""<<s_cwd<<"/FA_diagrams/subset_diagrams_"<<particle_name_tmp<<".pdf\",Paint[subdiags"<<particle_name_tmp<<"]];\n"  // print the FA diagram to pdf in local directory
  <<endl;
  }

  #ifdef RUN_ALL
  system("chmod +x make_figures.m ");
  if (verbose) system("./make_figures.m");
  else system("./make_figures.m  >/dev/null");
  #endif
  
  
  
}


// main routine to manage a diagram by diagram procedure


int main(int argc, char *argv[])
{

vector<std::string> particles, diagrams;
string model;

if (argc==1)
{
cout << "Please enter diagram number to calcuate and particle name, options are \"chi0\" and \"chi1\"  " << endl;
cout << "or alternatively enter the name to a list as \"-f list.txt\" " << endl;
}
else
{

string option = argv[1];
if (option=="-f")
{
  std::string particles_in [1000];
  std::string diagrams_in [1000]; int i=0;
  std::ifstream input(argv[2]);
  std::string line;
  while(getline(input, line))
  {
        if (!line.length() || line[0] == '#')
           continue;
        std::istringstream iss(line);
        if (i==0) iss >> model;
        else iss>> particles_in[i] >> diagrams_in[i];
        i=i+1;
  }
// run over all entries in the input file
cout << "drawing all diagrams specified in the input list" << endl;

  particles.resize(i);
  diagrams.resize(i);

  for (int k=0;k<i;k++)
  {
  particles[k]=particles_in[k];
  diagrams[k]=diagrams_in[k];
  }
  draw_diagrams(particles,diagrams,i,model);
  

}
else if (option == "-a")
{
//string particle;
//string model;
if (argc==4) {cout <<"drawing all diagrams for particle "<<argv[2]<<endl;;draw_all_diagrams(argv[2],argv[3]);}
else {cout << "please specify a particle after -a to use this option"<< endl;}


}


}


}
