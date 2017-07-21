/*
 Mass Builder

 James McKay
 Sep 2016
 Jan 2016

 --- utils.cpp ---

 This file contains functions used throughout the code to perform a range of
 tasks including reading from and writing to input/output streams, sorting
 lists of data and code generation

 */


#include "templates.hpp"
#include "cmake_variables.hpp"

namespace templates
{
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );

  void self_energy_hpp_header(ofstream &file)
  {
    file << "#ifndef SELF_ENERGY_H\n"
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
  }

  void data_input_block(ofstream &file)
  {

    file<<"  Data (){};\n"
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
  }


  void self_energy_src_preamble(ofstream &file, std::string particle_name, int subgroup)
  {
    file<< "/* ---------------------------------------------------- \n"
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
    << "#include \"cmake_variables.hpp\"\n"
    << "#include TSIL_PATH\n"
    << "  TSIL_DATA bar;\n"
    << " \n"
    << "#ifndef PI\n"
    << "#define PI 4.0L*atan(1.0L)\n"
    << "#endif\n"
    << "// define subroutines here\n"
    << "  TSIL_REAL Q2,Q;\n"
    << "  TSIL_REAL p;\n"
    << "  TSIL_COMPLEXCPP Log(TSIL_REAL a){complex<double> s(a/Q2,-0.000);return log(s);}\n"
    << "  TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}\n"
    << "  TSIL_COMPLEXCPP Power(TSIL_COMPLEXCPP a, int b){return pow(a,b);}\n"
    << "  TSIL_REAL Sin(TSIL_REAL a){return sin(a);}\n"
    << "  TSIL_REAL Cos(TSIL_REAL a){return cos(a);}\n"
    << "  TSIL_REAL Dot(TSIL_REAL a, TSIL_REAL b){return a*b;}\n"
    << "  TSIL_REAL Sqrt(TSIL_REAL a){return TSIL_POW(a,0.5);}\n"
    << "  TSIL_REAL Zeta2 = 1.6449340668482;\n"
    << "  TSIL_COMPLEXCPP Ae(TSIL_REAL a) { return TSIL_Aeps_(TSIL_POW(a,2),Q2);}\n"
    << "  TSIL_COMPLEXCPP Be(TSIL_REAL a, TSIL_REAL b) { return TSIL_Beps_(TSIL_POW(a,2),TSIL_POW(b,2), TSIL_POW(p,2), Q2);}\n"
    << "  int          init(Data data);\n"
    << "  TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}\n"
    << "  TSIL_COMPLEXCPP operator+(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}\n"
    << "  TSIL_COMPLEXCPP operator+(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c+b;}\n"
    << "  TSIL_COMPLEXCPP operator-(double a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}\n"
    << "  TSIL_COMPLEXCPP operator-(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c-b;}\n"
    << "  TSIL_COMPLEXCPP operator/(TSIL_COMPLEXCPP a,double b){TSIL_COMPLEXCPP c=b;return a/c;}\n"
    << "  TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}\n"
    << " \n";
  }

  void amp_preamble(ofstream &file)
  {
    file<< "/* ---------------------------------------------------- \n"
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
    << "#include \"cmake_variables.hpp\"\n"
    << "#include TSIL_PATH\n"
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
    << "\n";
  }


  void status_update(ofstream &file)
  {
    file << "Print[\" --------------------------------------- \"]\n"
    << "Print[\" The trial SE is:\"]\n"
    << "Print[\" --------------------------------------- \"]\n"
    << "Print[SEnTrial]\n"
    << "diff = Simplify[SEn-SEnTrial]\n"
    << "Print[\" --------------------------------------- \"]\n"
    << "Print[\" The difference between trial and actual SE is:\"]\n"
    << "Print[\" --------------------------------------- \"]\n"
    << "Print[diff]\n";
  }

  void user_input_guide()
  {
    cout << " Welcome to Mass Builder \n"
    << " ------------------------------------------------------------------------------------------------------------ \n"
    << " to use call ./mass_builder <flags> where required flags for each mode are given below:\n"
    << "                                                                                                      \n"
    << " -a -m <model>                                           -  compute all diagrams in models/<model>/diagrams.txt\n"
    << " -a -m <model>                               -i <file>   -  compute all diagrams in listed in file\n"
    << " -a -m <model>  -p <particle>  -d <diagram>              -  compute specific diagram\n"
    << " -g -m <model>                                           -  generate code for available diagrams\n"
    << " -g -m <model>                               -i <file>   -  generate code for diagrams listed in file\n"
    << " -f -m <model>  -p <particle>                            -  draw all FeynArts diagrams for particle\n"
    << " -e                                          -i <file>   -  evaluate self energy with values for masses and couplings in file\n"
    << " -x -m <model>                                           -  draw vertices listed in models/<model>/vertices\n"
    << " ------------------------------------- other additional flags ------------------------------------------------ \n"
    << " -z   -  use models/<model>/model.gen as the generic FeynArts model\n"
    << " -0   -  iteratively compute 1-loop pole masses for vector particles\n"
    << " -r   -  specify model restrictions, for example -r WinoLimit,WinoCouplings\n"
    << " ------------------------------------------------------------------------------------------------------------- " <<endl;
  }



  void print_math_header(std::string &input)
  {
    input += "SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True];";
    input += "SetOptions[DiracTrace, DiracTraceEvaluate -> True];";
    //input += "$GenericMixing = True;";
    input += "null=0;";
    input += "MassBuilderA[mass_, D_] := TAI[D, 0, {{1, mass}}];";
    input += "MassBuilderB[mass1_, mass2_, D_] := TBI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}}];";
    input += "MassBuilderJ[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}}];";
    input += "MassBuilderT[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{2, mass1}, {1, mass2}, {1, mass3}}];";
    input += "MassBuilderK[mass1_, mass2_, mass3_, D_] := TJI[D, 0, {{1, mass1}, {1, mass2}, {1, mass3}}];";
    input += "MassBuilderV[mass1_, mass2_, mass3_, mass4_, D_] := TVI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1, mass4}}];";
    input += "MassBuilderF[mass1_, mass2_, mass3_, mass4_, mass5_, D_] := TFI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1,mass4}, {1,mass5}}];";
  }



}
