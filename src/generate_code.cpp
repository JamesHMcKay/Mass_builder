/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 18/09/16
James McKay

--- generate_code.cpp ---

Creates the C++ interface to TSIL

make must be run in the build directory after this program has been run

*/

#include "utils.hpp"
#include "generate_code.hpp"

using namespace std;
using namespace utils;

namespace Generate_code
{


int get_loop_order(string type)
{
int loop_order;
string loop_order_str = utils::char_to_string(type[0]);
stringstream convert(loop_order_str);
convert >> loop_order;
return loop_order;
}

void generate_code (Options options)
{
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );
  
  
  // run a shell script to check if generator files exists, backup or create
  
  // we will take data from output/
  
  //  to start with we will create the functions.txt file
  
  
  string particle_name = ""; // hold this fixed, produce one cpp for each particle
  string particle_name_reduced = "";
  string tag = "empty";
  
  
  string coeff_integrals_tmp = "/output/coeff_integrals"; // vector containing file names
  const char* coeff_products_tmp = "/output/coeff_products"; // vector containing file names
  const char* summation_tmp = "/output/summation"; // vector containing file names
  
  
  
  //const char* coeff_integrals = coeff_integrals_tmp + tag + extn;
  
  
  const char *ext = ".txt";
  string underscore = "_";
  

  
  
  string model = options.model;
  
  
  
  
  ofstream functions;
  functions.open ("output/functions.txt");
  
  
  
  // need to get list of diagrams from a data file
  
  
  
  
  
  
  string c_file_diagrams = "";
  
  c_file_diagrams = options.input_list;
  
  const char *file_diagrams = c_file_diagrams.c_str();
  
  
  
  
  
  vector<std::string> tags;
  vector<std::string> particle_names,levels;
  string level;
  int nd; // number of diagrams, length of diagram_number vector
  
  cout << "input list = " << file_diagrams << endl;
  
  get_data(particle_names,tags,levels, nd,file_diagrams);
  
  
  
  
  
  for (int d = 0; d<nd;d++)
  {
  //if (particle_names[d] == particle_name_reduced)
  //{
  level = levels[d];
  tag = tags[d];
  particle_name = particle_names[d];
  
  
  particle_name =  part_name_simple(particle_names[d]);
  
  string c_coeff_integrals = "models/" + model +coeff_integrals_tmp + underscore + particle_name + underscore + tag + underscore + level + ext;
  string c_coeff_products = "models/" + model + coeff_products_tmp + underscore + particle_name + underscore + tag +underscore + level  + ext;
  string c_summation = "models/" + model + summation_tmp + underscore + particle_name + underscore + tag + underscore + level + ext;
  
  const char *coeff_integrals = c_coeff_integrals.c_str();
  const char *coeff_products = c_coeff_products.c_str();
  const char *summation = c_summation.c_str();
  
  
  
  
  
  functions<< "TSIL_COMPLEXCPP  diagram"<<"_"<< particle_name << "_" << tag << "_" << level << "()" <<endl;
  functions<< "{" << endl;
  
  
  
  ifstream infile(coeff_integrals);
  
  
  string content = "";
  int i;


  for(i=0 ; infile.eof()!=true ; i++) // get content of infile
  {
  content += infile.get();
  }
  i--;
  content.erase(content.end()-1);     // erase last character

 
  infile.close();

  functions << content;                 // output


  string content2 = "";

  ifstream infile2(coeff_products);

  for(i=0 ; infile2.eof()!=true ; i++) // get content of infile
  {
  content2 += infile2.get();
  }
  i--;
  content2.erase(content2.end()-1);     // erase last character

 
  infile2.close();

  functions << content2;                 // output

  string content3 = "";

  ifstream infile3(summation);

  for(i=0 ; infile3.eof()!=true ; i++) // get content of infile
  {
  content3 += infile3.get();
  }
  i--;
  content3.erase(content3.end()-1);     // erase last character

 
  infile3.close();

  functions << content3;                 // output
  functions << "}" << endl;
  functions << "     " << endl;
  functions << "     " << endl;
  functions << "     " << endl;


  /// end for loop over diagrams
  //}
  }
  
  

  functions.close();
  

  // -------------  now create the DoTSIL function, this should be done each time this is run to avoid duplicates -------------

  // need to know what what basis integrals are required, these are contained in the lists "integrals_tag.txt"


  int ni_total = 0;
  vector<std::string> integrals;
  
  for (int d = 0; d<nd;d++)
  {
  tag = tags[d];
  
  particle_name =  part_name_simple(particle_names[d]);
  
  level = levels[d];
  
 
  
  const char* file_integrals_tmp = "/output/basis_integrals"; // vector containing file names
  string c_file_integrals = "models/" + model + file_integrals_tmp + underscore + particle_name + underscore + tag + underscore + level + ext;
  const char *file_integrals = c_file_integrals.c_str();
  
  
  vector<std::string> integrals_temp;
  int ni; // number of diagrams, length of diagram_number vector
  
  
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
  
  string c_file_masses = "models/" + model+"/masses.txt";  // need to make this model independent
  const char *file_masses = c_file_masses.c_str();
  
  
  vector<string> masses_input,id_input;
  int na;
  get_data(masses_input,id_input,na,file_masses);
  std::map<std::string, Bases> base_map = set_bases(masses_input, id_input);
  
  
  
  
  
  
  
  
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
  << " \n"
  << "TSIL_DATA bar;\n"
  << " \n"
  << "#ifndef PI\n"
  << "#define PI 4.0L*atan(1.0L)\n"
  << "#endif\n"
  << "long double strtold(const char *, char **);\n"
  << "// define subroutines here\n"
  << "TSIL_REAL Power(TSIL_REAL a, int b){return TSIL_POW(a,b);}\n"
  << "TSIL_REAL Sin(TSIL_REAL a){return sin(a);}\n"
  << "TSIL_REAL Cos(TSIL_REAL a){return cos(a);}\n"
  << "TSIL_REAL Dot(TSIL_REAL a, TSIL_REAL b){return a*b;}\n"
  << "int          init(Data data);\n"
  << "TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}\n"
  << "TSIL_COMPLEXCPP Complex(double a,double b){dcomp i;i=-1;i=sqrt(i);TSIL_COMPLEXCPP result = a + i*b; return result ;}\n"
  <<endl;
  
  int count = 0;
  main_output<< "TSIL_COMPLEXCPP ";
  for (unsigned int i =0; i<integrals.size() ; i++)
  {
  count = count + 1;
  if (i!=(integrals.size()-1) && count !=5) {main_output << integrals[i] << ", ";}
  else{main_output << integrals[i] << " ";}
  if (count == 5) {count = 0; main_output << ";"<<endl; if (i!=(integrals.size()-1)){main_output<<"TSIL_COMPLEXCPP ";}}
  }
  if (count != 0){main_output<< ";"<<endl;}
  else {main_output<<endl;}
  
  
  // deal with Bxy = Byx case
  for (unsigned int i = 0; i<integrals.size();i++)
  {
  string name = integrals[i];
  
  Bases base_temp = base_map[name];
  base_temp.short_name = name;
  
  if (base_temp.type == "B")
  {
  Bases base_temp_B;
  base_temp_B.type = "B";
  base_temp_B.e1 = base_temp.e2;
  base_temp_B.e2 = base_temp.e1;
  string name_B = get_short_name(base_temp_B,masses_input, id_input);
  //integrals.push_back(name_B);
  main_output <<"TSIL_COMPLEXCPP " << name_B << ";" << endl;
  }
  }
  
  
  

  
  
  
  
  main_output << "TSIL_COMPLEXCPP  i;\n";
  vector<std::string> masses;
  int nm; // number of diagrams, length of diagram_number vector
  vector<std::string> temp_vec; // TODO remove!
  get_data(masses,temp_vec, nm,file_masses);
  main_output << "TSIL_REAL ";
  for (int i=0;i<nm;i++)
  {
  if (i!=(nm-1)) main_output << " " <<masses[i] << ", " << " " <<masses[i] << "2 , ";
  else main_output << " " <<masses[i] << ", " <<masses[i] << "2 ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
  // other global variable definitions
  
  
 main_output << "TSIL_REAL p,Pi;\n";
  // model specific variables need to be written here
  
  vector<std::string> couplings;
  string c_file_couplings = "models/" + model + "/couplings.txt";  // need to make this model independent
  const char *file_couplings = c_file_couplings.c_str();
  int nc; // number of diagrams, length of diagram_number vector
  get_data(couplings, nc,file_couplings);
  main_output << "TSIL_REAL ";
  for (int i=0;i<nc;i++)
  {
  if (i!=(nc-1)) main_output << couplings[i] << ", ";
  else main_output <<couplings[i] << " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";

  
  
  //  DoTSIL function
  
  
  main_output <<"void DoTSIL(TSIL_REAL s,TSIL_REAL Q2)\n"<<"{\n";
  

  
  for (unsigned int i = 0; i<integrals.size();i++)
  {
  string name = integrals[i];
  
  Bases base_temp = base_map[name];
  base_temp.short_name = name;
  
  print_doTSIL(main_output, base_temp);
  
  // deal with the Bxy = Byx case
  // this will be replaced by a more sophisticated algorithm which takes advantage of symmetries
  if (base_temp.type == "B")
  {
  Bases base_temp_B;
  base_temp_B.type = "B";
  base_temp_B.e1 = base_temp.e2;
  base_temp_B.e2 = base_temp.e1;
  string name_B = get_short_name(base_temp_B,masses_input, id_input);
  integrals.push_back(name_B);
  main_output << name_B << " = " << name << ";" << endl;
  
  }
  
  main_output << "\n";
  }
  
  
  
  main_output << "}"  << endl;
  main_output << "\n";
  main_output << "\n";
  // init function
  
  main_output << "  int init(Data data) \n"
  <<"{\n";

  for (int i=0;i<nc;i++)
  {
  if (i!=(nc-1)) main_output << couplings[i] << " = data."<<couplings[i]<<", ";
  else main_output << couplings[i] << " = data."<<couplings[i]<< " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";


  for (int i=0;i<nm;i++)
  {
  if (i!=(nm-1)) main_output <<" "<< masses[i] << " = data."<<masses[i]<<", " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) , ";
  else main_output << " "<<masses[i] << " = data."<<masses[i]<<" , " << masses[i]<<"2 = TSIL_POW(data."<<masses[i]<<", 2) ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";



  main_output<<"dcomp ii=-1;ii=sqrt(ii);i=ii;\n"
  <<"Pi=PI;\n"
  <<"return 0;\n"
  <<"}\n";
    
    
  
  
  
  
  ifstream infile4("output/functions.txt");
  
  string content2 = "";
  int ii;


  for(ii=0 ; infile4.eof()!=true ; ii++) // get content of infile
  {
  content2 += infile4.get();
  }
  ii--;
  content2.erase(content2.end()-1);     // erase last character

  
  infile4.close();

  main_output << content2;                 // output

  
  
  
  
  
  
  
  // MAIN FUNCTION CALL



  main_output << "void Self_energy::init_tsil(Data data)\n"
  << "{\n"
  << "TSIL_REAL qq= data.Q;\n"
  << "TSIL_REAL s=pow(data.P,2);\n"
  << "init(data);\n"
  << "DoTSIL(s,qq);\n"
  << "}\n"
  << "\n"
  << "\n"
  << "void Self_energy::run_tsil (Data &data) \n"
  << "{\n"
  << "p=data.P;\n"
  << "init_tsil(data);\n"
  << endl;

  vector<std::string> particle_names_short = particle_names;
  vector<std::string> particle_names_short_reduced;

  sort(particle_names_short.begin(),particle_names_short.end());
  particle_names_short.erase( unique( particle_names_short.begin(), particle_names_short.end() ), particle_names_short.end() );




  for (unsigned int i=0;i<particle_names_short.size();i++)
  {
    string particle_name_tmp = particle_names_short[i];
    string particle_name_tmp_short = part_name_simple(particle_names_short[i]);
    particle_names_short_reduced.push_back(particle_name_tmp_short);
    /* One loop */
    
    main_output<< "TSIL_COMPLEXCPP SE_1_"<<particle_name_tmp_short<<" = 0.L  ";
    for (int d = 0; d<nd;d++)
    {
    if (particle_names[d] == particle_name_tmp)
    {
    if (get_loop_order(levels[d]) == 1 )
    {
    
    main_output<< " + diagram"<<"_"<< particle_name_tmp_short << "_" << tags[d] << "_" << levels[d] << "()";
    
    tag = tags[d];
    }
    
    }
    }
    main_output << ";" << endl;
    main_output << "SE_1_"<<particle_name_tmp_short << " = " << "SE_1_"<<particle_name_tmp_short<<"*TSIL_POW(PI,2);"<<endl;
    //main_output << "cout << \"One-loop self energy of particle "<< particle_name_tmp_short << " = \" << real(SE_1_"<<particle_name_tmp_short<<") << endl;"<<endl;
    
    /* Two loop */
    main_output<< "TSIL_COMPLEXCPP SE_2_"<<particle_name_tmp_short<<" = 0.L ";
    for (int d = 0; d<nd;d++)
    {
    if (particle_names[d] == particle_name_tmp)
    {
    if (get_loop_order(levels[d]) == 2 )
    {
    
    main_output<< " + diagram"<<"_"<< particle_name_tmp_short << "_" << tags[d] << "_" << levels[d] << "()";
    
    tag = tags[d];
    }
    
    }
    }
    main_output << ";" << endl;
    main_output << "SE_2_"<<particle_name_tmp_short << " = " << "SE_2_"<<particle_name_tmp_short<<"*TSIL_POW(PI,4);"<<endl;
    //main_output << "cout << \"Two-loop self energy of particle "<< particle_name_tmp_short << " = \" << real(SE_2_"<<particle_name_tmp_short<<") << endl;"<<endl;
    
    
    //main_output << "data.SE_"<< particle_name_tmp_short << " = " << "real(SE_1_"<<particle_name_tmp_short<<" + " <<  "SE_2_"<<particle_name_tmp_short <<  ");"<<endl;
    main_output << "data.SE_1[\""<< particle_name_tmp_short << "\"] = " << "real(SE_1_"<<particle_name_tmp_short<< ");"<<endl;
    main_output << "data.SE_2[\""<< particle_name_tmp_short << "\"] = " << "real(SE_2_"<<particle_name_tmp_short << ");"<<endl;
  }

  main_output << "}" << endl;
  


  main_output.close();
  
  
  //////////////////////////////////////////////////////////////////////////////////////////
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

  //  first define all masses and couplings required by the self energy functions as doubles (they will later be recast to TSIL_REAL)


  for (int i=0;i<nc;i++)
  {
  data_h << "double "<< couplings[i]<<";"<<endl;
  }


  for (int i=0;i<nm;i++)
  {
  data_h << "double "<< masses[i]<<";"<<endl;
  }

  // create a variable for the self energy of each particle we have

 // for (unsigned int i=0;i<particle_names_short_reduced.size();i++)
 // {
 // data_h << "double SE_" << particle_names_short_reduced[i]<<";"<<endl;
 // }

  // create a map for self energies
  data_h << "std::map<std::string, double> SE_1;"<< endl;
  data_h << "std::map<std::string, double> SE_2;"<< endl;


  data_h << "double P, Q;" << endl;
  
  data_h << "std::vector<std::string> avail_part = {\"";
  for (unsigned int i=0;i<particle_names_short_reduced.size()-1;i++)
  {
  data_h << particle_names_short_reduced[i]<<"\",\"";
  }
  data_h << particle_names_short_reduced[particle_names_short_reduced.size()-1]<<"\"};"<<endl;

  
  

  // now create constructor and user input reader

  
    // constructor
  data_h<<"  Data (){};\n"
  <<"Data(Options options) {\n"
  <<"double param [99];\n"
  <<"std::string name [99]; int i=0;\n"
  <<"std::ifstream input(options.input_list);\n"
  <<"std::string line;\n"
  <<"while(getline(input, line)) {\n"
  <<"if (!line.length() || line[0] == '#')\n"
  <<" continue;\n"
  <<" std::istringstream iss(line);\n"
  <<" iss>> name[i] >> param[i];\n"
  <<"     i=i+1;\n"
  <<"   }\n"
  //<<"  }\n"
  <<"  for (int n=0;n<i+1;n++)\n"
  <<"  {\n";

  for (int i=0;i<nc;i++)
  {
  data_h<<"  if (name[n]==\""<<couplings[i] <<"\")\n"
  <<"  {\n"
  <<"  " <<couplings[i] << "=param[n];\n"
  <<"  }"<<endl;
  }

  for (int i=0;i<nm;i++)
  {
  data_h<<"  if (name[n]==\""<<masses[i] <<"\")\n"
  <<"  {"
  <<"  " <<masses[i] << "=param[n];"
  <<"  }"<<endl;
  }

  data_h<<"  if (name[n]==\"Q\")\n"
  <<"  {"
  <<"  Q =param[n];"
  <<"  }"<<endl;


  data_h<<"  if (name[n]==\"P\")\n"
  <<"  {"
  <<"  P =param[n];"
  <<"  }"<<endl;


  data_h<<"}\n"
  <<"}\n"
  <<"};\n"
  
  <<"#endif\n";

  
  
  
}
}
