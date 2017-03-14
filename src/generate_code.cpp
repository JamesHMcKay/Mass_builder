/*
 Mass Builder
 
 James McKay
 Aug - Sep 2016
 Jan 2017
 
 --- generate_code.cpp ---
 
 Creates the C++ interface to TSIL
 
 make must be run in the build directory after this program has been run
 */

#include "utils.hpp"
#include "generate_code.hpp"
#include "write_tsil_ini.hpp"

using namespace std;
using namespace utils;


void Generate_code::generate_self_energy_hpp()
{
  
  se_hpp<<"#ifndef SELF_ENERGY_H\n"
  <<"#define SELF_ENERGY_H\n"
  <<"#include \"data.hpp\"\n"
  <<"#include <complex>\n"
  <<"using namespace std;\n"
  <<"typedef complex<double> dcomp;\n"
  <<"class Self_energy\n"
  <<"{\n"
  <<"  private:\n"
  <<"  public:\n"
  <<"  Data data;\n"
  <<"  Self_energy (){}\n"
  <<"  void run_tsil(Data &data);\n"
  <<"};\n"
  <<"namespace tsil\n"
  <<"{\n"
  <<"  class Integrals\n"
  <<"  {\n"
  <<"    public:\n";
  
  decalare_var(se_hpp);
  
  se_hpp<<"    Integrals (){}\n"
  <<"    void DoTSIL(Data data);\n"
  <<"  };\n"
  <<"}\n";
}



void Generate_code::generate_data_hpp()
{
  
  
  // create header file with required masses and couplings which can be set via user input at runtime
  
  ofstream data_h;
  data_h.open ("include/data.hpp");
  
  
  data_h << "#ifndef DATA_H\n"
  <<"#define DATA_H\n"
  <<"#include \"options.hpp\"\n"
  <<"#include <map>\n"
  <<"using namespace std;\n"
  <<"struct Data\n"
  <<"{\n"
  <<"public:\n";
  
  // define all masses and couplings required by the self energy functions
  
  for (int i=0;i<nc;i++)
  {
    data_h << "  double "<< couplings[i]<<";"<<endl;
  }
  
  
  for (int i=0;i<nm;i++)
  {
    data_h << "  double "<< masses[i]<<";"<<endl;
  }
  
  
  // create a map for self energies
  data_h << "  std::map<std::string, double> SE_1;"<< endl;
  data_h << "  std::map<std::string, double> SE_2;"<< endl;
  data_h << "  double P, Q;" << endl;
  data_h << "  std::vector<std::string> avail_part = {\"";
  for (unsigned int i=0;i<particle_names_short_reduced.size()-1;i++)
  {
    data_h << particle_names_short_reduced[i]<<"\",\"";
  }
  data_h << particle_names_short_reduced[particle_names_short_reduced.size()-1]<<"\"};"<<endl;
  
  
  // create constructor and user input reader
  
  data_h<<"  Data (){};\n"
  <<"  Data(Options options) \n"
  <<"  {\n"
  <<"    double parameter [99];\n"
  <<"    std::string name [99]; int i=0;\n"
  <<"    std::ifstream input(options.input_list);\n"
  <<"    std::string line;\n"
  <<"    while(getline(input, line))\n"
  <<"    {\n"
  <<"      if (!line.length() || line[0] == '#')\n"
  <<"      continue;\n"
  <<"      std::istringstream iss(line);\n"
  <<"      iss>> name[i] >> parameter[i];\n"
  <<"      i=i+1;\n"
  <<"    }\n"
  <<"    for (int n=0;n<i+1;n++)\n"
  <<"    {\n";
  
  for (int i=0;i<nc;i++)
  {
    data_h<<"      if (name[n]==\""<<couplings[i] <<"\")\n"
    <<"      {\n"
    <<"        " <<couplings[i] << " = parameter[n];\n"
    <<"      }"<<endl;
  }
  
  for (int i=0;i<nm;i++)
  {
    data_h<<"      if (name[n]==\""<<masses[i] <<"\")\n"
    <<"      {\n"
    <<"        " <<masses[i] << " = parameter[n];\n"
    <<"      }"<<endl;
  }
  
  data_h<<"      if (name[n]==\"Q\")\n"
  <<"      {\n"
  <<"        Q = parameter[n];\n"
  <<"      }"<<endl;
  
  
  data_h<<"      if (name[n]==\"P\")\n"
  <<"      {\n"
  <<"        P = parameter[n];\n"
  <<"      }"<<endl;
  
  
  data_h<<"    }\n"
  <<"  }\n"
  <<"};\n"
  
  <<"#endif\n";
  data_h.close();
}


void Generate_code::decalare_var_tsil(ofstream &main_output)
{
  
  int count = 0;
  main_output<< "  TSIL_COMPLEXCPP ";
  for (unsigned int i =0; i<integrals.size() ; i++)
  {
    count = count + 1;
    if (i!=(integrals.size()-1) && count !=5) {main_output << integrals[i] << ", ";}
    else{main_output << integrals[i] << " ";}
    if (count == 5) {count = 0; main_output << ";"<<endl; if (i!=(integrals.size()-1)){main_output<<"  TSIL_COMPLEXCPP ";}}
  }
  if (count != 0){main_output<< ";"<<endl;}
  else {main_output<<endl;}
  
  // deal with Bxy = Byx case
  /*
   for (unsigned int i = 0; i<integrals.size();i++)
   {
   string name = integrals[i];
   
   Bases base_temp = base_map[name];
   base_temp.short_name = name;
   
   if (base_temp.type == "B" && base_temp.e1!=base_temp.e2)
   {
   Bases base_temp_B;
   base_temp_B.type = "B";
   base_temp_B.e1 = base_temp.e2;
   base_temp_B.e2 = base_temp.e1;
   string name_B = get_short_name(base_temp_B,masses_input, id_input);
   main_output <<"  TSIL_COMPLEXCPP " << name_B << ";" << endl;
   }
   }
   */
  
  main_output << "  TSIL_COMPLEXCPP  i;\n";
  
  
  
  get_data(masses,temp_vec, nm,file_masses);
  main_output << "  TSIL_REAL ";
  for (int i=0;i<nm;i++)
  {
    if (i!=(nm-1)) main_output << " " <<masses[i] << ", " << " " <<masses[i] << "2 , ";
    else main_output << " " <<masses[i] << ", " <<masses[i] << "2 ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
  // other global variable definitions
  
  
  main_output << "  TSIL_REAL p,Pi;\n";
  // model specific variables need to be written here
  
  
  string c_file_couplings = "models/" + model + "/couplings.txt";  // need to make this model independent
  const char *file_couplings = c_file_couplings.c_str();
  
  get_data(couplings, relationships, nc,file_couplings,true);
  main_output << "  TSIL_REAL ";
  
  for (int i=0;i<nc;i++)
  {
    if (i!=(nc-1)) main_output << couplings[i] << ", ";
    else main_output <<couplings[i] << " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
}

void Generate_code::decalare_var(ofstream &main_output)
{
  
  int count = 0;
  main_output<< "    std::complex<long double> ";
  for (unsigned int i =0; i<integrals.size() ; i++)
  {
    count = count + 1;
    if (i!=(integrals.size()-1) && count !=5) {main_output << integrals[i] << ", ";}
    else{main_output << integrals[i] << " ";}
    if (count == 5) {count = 0; main_output << ";"<<endl; if (i!=(integrals.size()-1)){main_output<<"    std::complex<long double> ";}}
  }
  if (count != 0){main_output<< ";"<<endl;}
  else {main_output<<endl;}
  
  // deal with Bxy = Byx case
  for (unsigned int i = 0; i<integrals.size();i++)
  {
    string name = integrals[i];
    
    Bases base_temp = base_map[name];
    base_temp.short_name = name;
    
    if (base_temp.type == "B" && base_temp.e1!=base_temp.e2)
    {
      Bases base_temp_B;
      base_temp_B.type = "B";
      base_temp_B.e1 = base_temp.e2;
      base_temp_B.e2 = base_temp.e1;
      string name_B = get_short_name(base_temp_B,masses_input, id_input);
      main_output <<"    std::complex<long double> " << name_B << ";" << endl;
    }
  }
  
  main_output << "    std::complex<long double>  i;\n";
  
  
  
  get_data(masses,temp_vec, nm,file_masses);
  main_output << "    double ";
  for (int i=0;i<nm;i++)
  {
    if (i!=(nm-1)) main_output << " " <<masses[i] << ", " << " " <<masses[i] << "2 , ";
    else main_output << " " <<masses[i] << ", " <<masses[i] << "2 ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
  // other global variable definitions
  
  
  main_output << "    double p,Pi;\n";
  // model specific variables need to be written here
  
  
  string c_file_couplings = "models/" + model + "/couplings.txt";  // need to make this model independent
  const char *file_couplings = c_file_couplings.c_str();
  
  get_data(couplings, relationships, nc,file_couplings,true);
  main_output << "  double ";
  
  for (int i=0;i<nc;i++)
  {
    if (i!=(nc-1)) main_output << couplings[i] << ", ";
    else main_output <<couplings[i] << " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
}


void Generate_code::generate_particle_src(std::string particle,int subgroup)
{
  particle_name =  part_name_simple(particle);
  string c_src_file_name = "src/amp_" + particle_name + "_" + std::to_string(subgroup) + cpp;
  
  const char *src_file_name = c_src_file_name.c_str();
  
  ofstream functions;
  functions.open (src_file_name);
  
  cout << "writing source file = " << src_file_name << endl;
  
  
  /* TSIL_PATH */ std::string TSIL = "/Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h";
  
  functions<< "/* ---------------------------------------------------- \n"
  << "(* This file has been automatically generated by generate_code.cpp, on "<< now->tm_mday << '-'
  << (now->tm_mon + 1) << '-'<< (now->tm_year + 1900) <<", do not edit *)\n"
  << " ---------------------------------------------------- */ \n"
  << "#include <iostream>\n"
  << "#include <iomanip>\n"
  << "#include <string>\n"
  << "#include <fstream>\n"
  << "#include <sstream>\n"
  << "#include <complex>\n"
  << "#include <vector>\n"
  << "#include <cmath>\n"
  << "#include <cfloat>\n"
  << "#include <cstdlib>\n"
  << "#include <ctime>\n"
  << "#include \"self_energy.hpp\"\n"
  << " \n"
  << "using namespace std;\n"
  << "using namespace tsil;\n"
  << " \n"
  << " \n"
  << " \n"
  << "namespace " << particle_name << "_" << std::to_string(subgroup) << "\n"
  << "{\n"
  << "#include \""<< TSIL << "\"\n"
  << "  TSIL_DATA bar;\n"
  << " \n"
  << "#ifndef PI\n"
  << "#define PI 4.0L*atan(1.0L)\n"
  << "#endif\n"
  << "// define subroutines here\n"
  << "  TSIL_REAL Q2,Q;\n"
  << "  TSIL_REAL Log(TSIL_REAL a){return TSIL_LOG(a/Q2);}\n"
  << "  TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}\n"
  << "  TSIL_REAL Sin(TSIL_REAL a){return sin(a);}\n"
  << "  TSIL_REAL Cos(TSIL_REAL a){return cos(a);}\n"
  << "  TSIL_REAL Dot(TSIL_REAL a, TSIL_REAL b){return a*b;}\n"
  << "  TSIL_REAL Zeta2 = 1.6449340668482;\n"
  << "  TSIL_COMPLEXCPP Ae(TSIL_REAL a) { return TSIL_Aeps_(a,Q2);}\n"
  << "  TSIL_COMPLEXCPP Be(TSIL_REAL a, TSIL_REAL b, TSIL_REAL p) { return TSIL_Beps_(a,b, TSIL_POW(p,2), Q2);}\n"
  << "  int          init(Data data);\n"
  << "  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}\n"
  << "  TSIL_COMPLEXCPP operator+(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}\n"
  << "  TSIL_COMPLEXCPP operator-(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}\n"
  << "  TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}\n"
  << " \n"
  <<endl;
  
  
  decalare_var_tsil(functions);
  
  functions << "\n";
  functions << "\n";
  functions << "\n";
  
  for (int d = 0; d<nd;d++)
  {
    if (particle == particle_names[d] && subgroup==subgrouplist[d])
    {
      
      level = levels[d];
      tag = tags[d];
      
      particle_name =  part_name_simple(particle_names[d]);
      
      string c_coeff_integrals = "models/" + model +coeff_integrals_tmp  + particle_name + underscore + tag + underscore + level + ext;
      string c_coeff_products = "models/" + model + coeff_products_tmp  + particle_name + underscore + tag +underscore + level  + ext;
      string c_summation = "models/" + model + summation_tmp + particle_name + underscore + tag + underscore + level + ext;
      
      const char *coeff_integrals = c_coeff_integrals.c_str();
      const char *coeff_products = c_coeff_products.c_str();
      const char *summation = c_summation.c_str();
      
      functions<< "TSIL_COMPLEXCPP  diagram"<<"_"<< particle_name << "_" << tag << "_" << level << "()" <<endl;
      functions<< "{" << endl;
      
      ifstream infile(coeff_integrals);
      
      string content = "";
      int i;
      
      for(i=0 ; infile.eof()!=true ; i++)
      {
        content += infile.get();
      }
      i--;
      content.erase(content.end()-1);
      
      infile.close();
      
      functions << content;
      
      string content2 = "";
      
      ifstream infile2(coeff_products);
      
      for(i=0 ; infile2.eof()!=true ; i++)
      {
        content2 += infile2.get();
      }
      i--;
      content2.erase(content2.end()-1);
      
      infile2.close();
      
      functions << content2;
      
      string content3 = "";
      
      ifstream infile3(summation);
      
      for(i=0 ; infile3.eof()!=true ; i++)
      {
        content3 += infile3.get();
      }
      i--;
      content3.erase(content3.end()-1);
      
      
      infile3.close();
      
      functions << content3;
      functions << "}" << endl;
      functions << "     " << endl;
      functions << "     " << endl;
      functions << "     " << endl;
    }
    
    
    
    
  }
  
  
  functions << "  void SE_"<<particle_name << " (Data data, Integrals integral)\n";
  functions << "  {\n";
  
  for (int i=0;i<nc;i++)
  {
    if (i!=(nc-1)) functions << "    " << couplings[i] << " = data."<<couplings[i]<<", ";
    else functions << "    " << couplings[i] << " = data."<<couplings[i]<< " ";
  }
  functions<<";\n"
  << "\n"
  << "    Q2 = data.Q;\n"
  << "    Q = pow(Q2,0.5);\n"
  << "    p=data.P;\n";
  
  
  for (int i=0;i<nm;i++)
  {
    if (i!=(nm-1)) functions <<"    "<< masses[i] << " = data."<<masses[i]<<", " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) , ";
    else functions << "    "<<masses[i] << " = data."<<masses[i]<<" , " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) ";
  }
  functions<<";\n"
  << "\n";
  
  
  
  functions<<"    dcomp ii=-1;ii=sqrt(ii);i=ii;\n"
  <<"    Pi=PI;\n";
  
  if (relationships.size()!=0)
  {
    for (unsigned int i=0;i<relationships.size();i++)
    {
      functions << "    " << couplings[i] << " = " << relationships[i] << ";"<<endl;
    }
  }
  
  
  // set values of integrals
  
  for (unsigned int i =0; i<integrals.size() ; i++)
  {
    functions << "    " << integrals[i] << " = integral."<<integrals[i]<<";\n";
  }
  
  // deal with Bxy = Byx case
  for (unsigned int i = 0; i<integrals.size();i++)
  {
    string name = integrals[i];
    Bases base_temp = base_map[name];
    base_temp.short_name = name;
    if (base_temp.type == "B" && base_temp.e1!=base_temp.e2)
    {
      Bases base_temp_B;
      base_temp_B.type = "B";
      base_temp_B.e1 = base_temp.e2;
      base_temp_B.e2 = base_temp.e1;
      string name_B = get_short_name(base_temp_B,masses_input, id_input);
      functions <<  "    " << name_B << " = integral." << name_B << ";\n";
    }
  }
  // done setting integral values
  
  functions << "  }\n"
  << "\n"
  << "\n"
  << "  double SE_1()\n"
  << "  {\n";
  
  /* One loop */
  
  functions<< "    double SE = real(0.L  ";
  for (int d = 0; d<nd;d++)
  {
    if (particle_names[d] == particle && subgroup==subgrouplist[d])
    {
      if (get_loop_order(levels[d]) == 1 )
      {
        
        functions<< " + diagram"<<"_"<< particle_name << "_" << tags[d] << "_" << levels[d] << "()";
        
        tag = tags[d];
      }
      
    }
  }
  functions <<");\n"
  << "    SE  = -SE*TSIL_POW(PI,2);\n"
  << "    return SE;\n"
  <<"  }\n"
  << "\n"
  << "\n"
  << "  double SE_2()\n"
  << "  {\n";
  
  /* Two loop */
  functions<< "    double  SE = real(0.L ";
  for (int d = 0; d<nd;d++)
  {
    if (particle_names[d] == particle && subgroup==subgrouplist[d])
    {
      if (get_loop_order(levels[d]) == 2 )
      {
        
        functions<< " + diagram"<<"_"<< particle_name << "_" << tags[d] << "_" << levels[d] << "()";
        
        tag = tags[d];
      }
      
    }
  }
  functions << ");" << endl;
  functions << "    SE = -SE * TSIL_POW(PI,4);\n"<<endl;
  functions << "    return SE;\n"<<endl;
  functions <<"  }\n";
  functions << "\n"
  <<"}\n";
  
  functions.close();
  
  se_hpp << "namespace " << particle_name << "_" << std::to_string(subgroup) << "\n"
  <<"{\n"
  <<" \n"
  <<"void  SE_"<<particle_name <<"(Data data, tsil::Integrals integral);\n"
  <<"\n"
  <<"double SE_1();"
  <<"double SE_2();"
  <<"\n"
  <<"}\n";
  //    Print out value of each amplitude to terminal  (REDUNDANT OPTION NOW, UPDATE LATER)
  //    if (options.detailed_output)
  //    {
  //      for (int d = 0; d<nd;d++)
  //      {
  //        if (particle_names[d] == particle_name_tmp)
  //        {
  //          if (get_loop_order(levels[d]) == 2 )
  //          {
  //            main_output<< "  cout << \"diagram"<<" "<< particle_name_tmp_short << "_" <<  tags[d] << " = \" << TSIL_POW(PI,4)*real(diagram" <<"_"<< particle_name_tmp_short << "_" << tags[d] << "_" << levels[d] << "())<<endl;"<<endl;
  //          }
  //        }
  //      }
  //    }
}



void Generate_code::generate_code()
{
  
  // sorting integrals
  int ni_total = 0;
  
  for (int d = 0; d<nd;d++)
  {
    tag = tags[d];
    particle_name =  part_name_simple(particle_names[d]);
    level = levels[d];
    const char* file_integrals_tmp = "/output/basis_integrals"; // vector containing file names
    string c_file_integrals = "models/" + model + file_integrals_tmp + underscore + particle_name + underscore + tag + underscore + level + ext;
    const char *file_integrals = c_file_integrals.c_str();
    vector<std::string> integrals_temp;
    int ni; // number of diagrams
    get_data(integrals_temp, ni,file_integrals);
    integrals.resize(ni_total+ni);
    for (int i =0; i<ni ; i++)
    {
      integrals[ni_total+i] = integrals_temp[i];
    }
    ni_total = ni_total + ni;
  }
  // we have a list of integrals so now need to remove duplicates and pass to routines for determining the required TSIL functions
  sort(integrals.begin(),integrals.end());
  integrals.erase( unique( integrals.begin(), integrals.end() ), integrals.end() );
  string c_file_masses = "models/" + model+"/masses.txt";
  file_masses = c_file_masses.c_str();
  get_data(masses_input,id_input,na,file_masses);
  base_map = set_bases(masses_input, id_input);
  
  
  generate_self_energy_hpp();
  
  // MAIN HEADER OF FILE
  
  
  //const char *ext_cpp = ".cpp";
  const char* main_output_tmp = "src/self_energy.cpp"; // vector containing file names
  string c_main_output = main_output_tmp;// + ext_cpp;
  const char *main_output_file = c_main_output.c_str();
  
  ofstream main_output;
  main_output.open (main_output_file);
  
  /* TSIL_PATH */ std::string TSIL = "/Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h";
  
  main_output<< "/* ---------------------------------------------------- \n"
  << "(* This file has been automatically generated by generate_code.cpp, on "<< now->tm_mday << '-'
  << (now->tm_mon + 1) << '-'<< (now->tm_year + 1900) <<", do not edit *)\n"
  << " ---------------------------------------------------- */ \n"
  << "#include <iostream>\n"
  << "#include <iomanip>\n"
  << "#include <string>\n"
  << "#include <fstream>\n"
  << "#include <sstream>\n"
  << "#include <complex>\n"
  << "#include <vector>\n"
  << "#include <cmath>\n"
  << "#include <cfloat>\n"
  << "#include <cstdlib>\n"
  << "#include <ctime>\n"
  << "#include \""<< TSIL << "\"\n"
  << "#include \"self_energy.hpp\"\n"
  << " \n"
  << "using namespace std;\n"
  << " \n"
  << " \n"
  << "namespace tsil"
  << "{\n"
  << "  TSIL_DATA bar;\n"
  << " \n"
  << "#ifndef PI\n"
  << "#define PI 4.0L*atan(1.0L)\n"
  << "#endif\n"
  <<endl;
  
  
  
  
  //  DoTSIL function
  
  
  
  main_output <<" void Integrals::DoTSIL(Data data)\n"<<"{\n";
  // init function
  
  for (int i=0;i<nc;i++)
  {
    if (i!=(nc-1)) main_output << "  " << couplings[i] << " = data."<<couplings[i]<<", ";
    else main_output << "  " << couplings[i] << " = data."<<couplings[i]<< " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  main_output<< "   TSIL_REAL Q2 = data.Q;\n";
  main_output<< "   TSIL_REAL Q = pow(Q2,0.5);\n";
  main_output<< "   TSIL_REAL s=pow(data.P,2);\n";
  
  
  for (int i=0;i<nm;i++)
  {
    if (i!=(nm-1)) main_output <<"  "<< masses[i] << " = data."<<masses[i]<<", " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) , ";
    else main_output << "  "<<masses[i] << " = data."<<masses[i]<<" , " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
  
  
  main_output<<"    dcomp ii=-1;ii=sqrt(ii);i=ii;\n"
  <<"   Pi=PI;\n";
  
  if (options.optimise)
  {
    Print_dotsil print_tsil(integrals, base_map);
    print_tsil.print_to_file(main_output);
  }
  else
  {
    int ni = integrals.size();
    for (int i = 0; i<ni;i++)
    {
      string name = integrals[i];
      
      Bases base_temp = base_map[name];
      base_temp.short_name = name;
      
      print_doTSIL(main_output, base_temp);
      
      // deal with the Bxy = Byx case
      if (base_temp.type == "B" && base_temp.e1!=base_temp.e2)
      {
        Bases base_temp_B;
        base_temp_B.type = "B";
        base_temp_B.e1 = base_temp.e2;
        base_temp_B.e2 = base_temp.e1;
        string name_B = get_short_name(base_temp_B,masses_input, id_input);
        integrals.push_back(name_B);
        main_output << "  " <<  name_B << " = " << name << ";" << endl;
      }
      main_output << "\n";
    }
    
  }
  
  main_output << "  }"  << endl;
  main_output << "}"  << endl;
  main_output << "\n";
  main_output << "\n";
  
  
  main_output << "void Self_energy::run_tsil (Data &data)\n"
  <<"{\n"
  <<"  tsil::Integrals integrals;\n"
  <<"  integrals.DoTSIL(data);\n";
  
  // sort list into groups
  // for each particle, call function which will generate the source file
  // can soon split this into smaller subgroups
  
  vector<std::string> particle_names_short = particle_names;
  
  
  sort(particle_names_short.begin(),particle_names_short.end());
  particle_names_short.erase( unique( particle_names_short.begin(), particle_names_short.end() ), particle_names_short.end() );
  
  
  // for each particle split the list of diagrams into subgroups of 10 and create a new particle tag P1_1, P1_2, P1_3, ...
  
  subgrouplist.resize(tags.size());
  
  for (unsigned int i=0;i<particle_names_short.size();i++)
  {
    string particle_name_tmp = particle_names_short[i];
    string particle_name_tmp_short = part_name_simple(particle_names_short[i]);
    particle_names_short_reduced.push_back(particle_name_tmp_short);
    
    
    // function to determine number of diagrams associated with this particle -> nd (number of groups of 10, remainder)
    pair<int,int> num_diagrams = number_of_diagrams(particle_names_reduced[i]);
    
    
    main_output << "  double "<<  particle_name_tmp_short << "_1= 0;\n";
    main_output << "  double "<<  particle_name_tmp_short << "_2= 0;\n";
    
    // loop over nd/10 setting setting a flag for include diagram or not (new vector of bool)
    
    for (int j = 0; j< num_diagrams.first+1; j++)
    {
      generate_particle_src(particle_names_reduced[i],j);
      main_output<< "  " << particle_name_tmp_short <<"_"<<std::to_string(j)<<"::SE_"<<particle_name_tmp_short<<"(data,integrals);\n";
      main_output << "  "<<  particle_name_tmp_short << "_1 = " << particle_name_tmp_short << "_1 + " << particle_name_tmp_short <<"_"<<std::to_string(j)<<"::SE_1();\n";
      main_output << "  "<<  particle_name_tmp_short << "_2 = " << particle_name_tmp_short << "_2 + " << particle_name_tmp_short <<"_"<<std::to_string(j)<<"::SE_2();\n";
    }
  
    
    //main_output << "  double "<<  particle_name_tmp_short << "_1= " << particle_name_tmp_short <<"::SE_1();\n";
    //main_output << "  double "<<  particle_name_tmp_short << "_2= " << particle_name_tmp_short <<"::SE_2();\n";
    // this step should call all available functions tagged with this particle name and loop level
    
    main_output << "  data.SE_1[\""<< particle_name_tmp_short << "\"] = " << "real("<<particle_name_tmp_short<< "_1);"<<endl;
    main_output << "  data.SE_2[\""<< particle_name_tmp_short << "\"] = " << "real("<<particle_name_tmp_short << "_2);"<<endl;
  }
  main_output<<"}\n";
  main_output.close();
  
  
  
  generate_data_hpp();
  
  se_hpp<<"#endif\n";
  se_hpp.close();
  
  
  
  
  
  
}
