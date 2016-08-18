/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 18/09/16
James McKay

--- generate_code.cpp ---

Creates the C++ interface to TSIL using the data in Mass_builder/Generator/output/,
the resultant file is put in Mass_builder/src/self_energy.cpp.

make must be run in the build directory after this program has been run

*/




/* goal: generate cpp code in the directory ../src/ that contains a function
definition containing the coefficients and the total self energy SE = ...
also create a DoTSIL function to evaluate the integrals

function should have generic name with suffix given by particle name
and diagram number

need capacity to add multiple such functions to the file and update (removing
any duplications) the DoTSIL function.

---------

--- the process ---

need three files:

 - integrals.txt
 - functions.txt
 - diagrams.txt

these will be stored in output/generator/
should give option to start fresh (but in doing so previous versions
should be backed up to output/generator/backup/)

for a new diagram we do the following:

1.) add the required integrals to integrals.txt -- sort and remove duplicates
2.) add the new functions, containing all the coefficients to functions.txt
3.) add the new diagram number to diagrams.txt

4.) generate new .cpp code, create fresh DoTSIL function, straight copy of 
functions.txt and then a fresh master call, which will call all diagrams
listed in diagrams.txt


*/


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



void get_data_2_columns(vector<std::string> &A, vector<std::string> &B,int &n,const char *filename)
{

cout << "reading file = " << filename << endl;

std::ifstream input(filename);
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
A.resize(n);
B.resize(n);


n=0;
std::ifstream input2(filename);
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  iss2>> A[n] >> B[n];
    n=n+1;
 }
input.close();
input2.close();
}

void get_data(vector<std::string> &A,int &n,const char *filename)
{

cout << "reading file = " << filename << endl;

std::ifstream input(filename);
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
A.resize(n);


n=0;
std::ifstream input2(filename);
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  iss2>> A[n];
    n=n+1;
 }
 input.close();
 input2.close();
}







void print_TSIL_A(ofstream &myfile, string elements)
{
myfile << "A"<<elements << " = -i*TSIL_A_ (m"<<elements<<" , Q2);"<<endl;
}

void print_TSIL_B(ofstream &myfile, string elements)
{

myfile << "B"<< elements <<" = i*TSIL_B_ (m" << elements[0] << ", m" << elements[1] << ", s, Q2);"<< endl;
}

void print_TSIL_K(ofstream &myfile, string elements)
{

myfile << "K"<< elements <<" = TSIL_I2_(m" << elements[0] << ", m" << elements[1] << ", m" << elements[2] << ", Q2);"<< endl;
}


void print_TSIL_J(ofstream &myfile, string elements)
{
  char A = elements[0],B=elements[1],C=elements[2];
  //string output0 = A+B+C;
  
  
  myfile << "TSIL_SetParametersST (&bar,m" << B << ", m" << A << ", m" << C <<", Q2);" << endl;  // write TSIL evaluate statement for this name sequence
  myfile << "TSIL_Evaluate (&bar, s);" << endl;
  myfile << "J" << elements << "= TSIL_GetFunction (&bar,\"Suxv" <<"\");"<< endl;    // TSIL get function statement
  
}


void print_TSIL_T(ofstream &myfile, string elements)
{
  char A = elements[0],B=elements[1],C=elements[2];
  //string output0 = A+B+C;
  
  
  myfile << "TSIL_SetParametersST (&bar,m" << A << ", m" << B << ", m" << C <<", Q2);" << endl;  // write TSIL evaluate statement for this name sequence
  myfile << "TSIL_Evaluate (&bar, s);" << endl;
  myfile << "T" << elements << "= -TSIL_GetFunction (&bar,\"Txuv" <<"\");"<< endl;    // TSIL get function statement
  
}


void print_TSIL_F(ofstream &myfile, string elements)
{
  char A = elements[0],B=elements[1],C=elements[2],D=elements[3],E=elements[4];
  //string output0 = A+B+C+D+E;

  myfile << "TSIL_SetParameters (&bar,m" << A << ", m" << B << ", m" << C << " , m" << D << " , m" << E  << ", Q2);" << endl;  // write TSIL evaluate statement for this name sequence
  myfile << "TSIL_Evaluate (&bar, s);" << endl;
  myfile << "F" << elements << "= TSIL_GetFunction (&bar,\"M" <<"\");"<< endl;    // TSIL get function statement

}





void print_TSIL_V(ofstream &myfile, string elements)
{

  char A = elements[0],B=elements[1],C=elements[2],D=elements[3];

  myfile << "TSIL_SetParameters (&bar,m" << D << ", m" << C << ", m" << B << " , " << "1.0" << " , m" << A  << ", Q2);" << endl;  // write TSIL evaluate statement for this name sequence
  myfile << "TSIL_Evaluate (&bar, s);" << endl;
  myfile << "V" << elements << "= -TSIL_GetFunction (&bar,\"Uzxyv" <<"\");"<< endl;    // TSIL get function statement

}











int main (int argc, char *argv[])
{
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );
  
  
  // run a shell script to check if generator files exists, backup or create
  
  // we will take data from output/
  
  //  to start with we will create the functions.txt file
  
  
  string particle_name = "_chi0_"; // hold this fixed, produce one cpp for each particle
  string particle_name_reduced = "chi0";
  //int diagram_number = 13; // make this into a vector of diagrams
  string tag = "empty";
  
  
  string coeff_integrals_tmp = "output/coeff_integrals"; // vector containing file names
  const char* coeff_products_tmp = "output/coeff_products"; // vector containing file names
  const char* summation_tmp = "output/summation"; // vector containing file names
  
  
  
  //const char* coeff_integrals = coeff_integrals_tmp + tag + extn;
  
  
  const char *ext = ".txt";
  string underscore = "_";
  

  
  
  
  
  
  
  ofstream functions;
  functions.open ("generator/functions.txt");
  
  
  
  // need to get list of diagrams from a data file
  
  
  
  
  string c_file_diagrams = "output/diagrams.txt";
  if (argc==1)
  {
  
  cout << "using all available diagrams"<< endl;
  }
  else {
  
  c_file_diagrams = argv[1];
  
  }
  
  
  
  
  
  const char *file_diagrams = c_file_diagrams.c_str();
  
  
  
  
  
  vector<std::string> tags;
  vector<std::string> particle_names;
  int nd; // number of diagrams, length of diagram_number vector
  get_data_2_columns(particle_names,tags, nd,file_diagrams);
  
  
  
  for (int d = 0; d<nd;d++)
  {
  //if (particle_names[d] == particle_name_reduced)
  //{
  
  tag = tags[d];
  particle_name = particle_names[d];
  
  string c_coeff_integrals = coeff_integrals_tmp + underscore + particle_name + underscore + tag + ext;
  string c_coeff_products = coeff_products_tmp + underscore + particle_name + underscore + tag + ext;
  string c_summation = summation_tmp + underscore + particle_name + underscore + tag + ext;
  
  const char *coeff_integrals = c_coeff_integrals.c_str();
  const char *coeff_products = c_coeff_products.c_str();
  const char *summation = c_summation.c_str();
  
  
  
  
  
  functions<< "TSIL_COMPLEXCPP  diagram"<<"_"<< particle_name << "_" << tag << "()" <<endl;
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

  cout << i << " characters read...\n";
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

  cout << i << " characters read...\n";
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

  cout << i << " characters read...\n";
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
  //if (particle_names[d] == particle_name_reduced)
 // {
  
  const char* file_integrals_tmp = "output/basis_integrals"; // vector containing file names
  string c_file_integrals = file_integrals_tmp + underscore + particle_names[d] + underscore + tag + ext;
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
  //}

  }
  
  
  // we have a list of integrals so now need to remove duplicates and pass to routines for determining the require TSIL functions
  
  sort(integrals.begin(),integrals.end());
  integrals.erase( unique( integrals.begin(), integrals.end() ), integrals.end() );
  
  
  
  
  
  // MAIN HEADER OF FILE

  
  //const char *ext_cpp = ".cpp";
  const char* main_output_tmp = "../src/self_energy.cpp"; // vector containing file names
  string c_main_output = main_output_tmp;// + ext_cpp;
  const char *main_output_file = c_main_output.c_str();
  
  ofstream main_output;
  main_output.open (main_output_file);
  
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
  << "#include \"/Users/jamesmckay/Documents/Programs/tsil-1.3/tsil_cpp.h\"  // Required TSIL header file\n"
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
  << "int          init(Data data);\n"
  << "TSIL_COMPLEXCPP operator*(int a, TSIL_COMPLEXCPP b){TSIL_COMPLEXCPP c=a;return c*b;}\n"
  <<endl;
  
  int count = 0;
  main_output<< "TSIL_COMPLEXCPP ";
  for (int i =0; i<integrals.size() ; i++)
  {
  count = count + 1;
  if (i!=(integrals.size()-1) && count !=5) {main_output << integrals[i] << ", ";}
  else{main_output << integrals[i] << " ";}
  if (count == 5) {count = 0; main_output << ";"<<endl; if (i!=(integrals.size()-1)){main_output<<"TSIL_COMPLEXCPP ";}}
  }
  if (count != 0){main_output<< ";"<<endl;}
  else {main_output<<endl;}
  
  main_output << "TSIL_COMPLEXCPP  i=Power(-1,0.5);\n";
  vector<std::string> masses;
  string c_file_masses = "list_input.txt";
  const char *file_masses = c_file_masses.c_str();
  int nm; // number of diagrams, length of diagram_number vector
  get_data(masses, nm,file_masses);
  main_output << "TSIL_REAL ";
  for (int i=0;i<nm;i++)
  {
  if (i!=(nm-1)) main_output << "m" <<masses[i] << ", ";
  else main_output << "m" <<masses[i] << " ";
  }
  main_output<<";"<<endl;
  main_output<< "\n";
  
  // other global variable definitions
 main_output << "TSIL_REAL Pi,g2, tW, sw2, cw2, sw, cw , C,p,S2TW;\n";
  
  
  
  
  
  //  DoTSIL function
  
  
  main_output <<"void DoTSIL(TSIL_REAL s,TSIL_REAL Q2)\n"<<"{\n";
  
  // sort the integrals into types and pass each list of each type to a seperate sorting and printing routine
  
  //vector<std::string> V_A,V_B,V_C,V_D;
  //int nv = 0;
  
  for (int i = 0; i<integrals.size();i++)
  {
  string name = integrals[i];
  stringstream _e1,_e2,_e3,_e4,_e5, _type;
  string e1,e2,e3,e4,e5, type;
  _type << name[0];
  _type >> type;
  _e1 << name[1];_e1 >> e1;
  _e2 << name[2];_e2 >> e2;
  _e3 << name[3];_e3 >> e3;
  _e4 << name[4];_e4 >> e4;
  _e5 << name[5];_e5 >> e5;
  
  if (type == "A") print_TSIL_A(main_output,e1);
  if (type == "B") print_TSIL_B(main_output, e1+e2);
  if (type == "J") print_TSIL_J(main_output, e1+e2+e3);
  if (type == "T") print_TSIL_T(main_output, e1+e2+e3);
  if (type == "K") print_TSIL_K(main_output, e1+e2+e3);
  if (type == "V") print_TSIL_V(main_output,e1+e2+e3+e4);//{V_A.push_back(e1);V_B.push_back(e2);V_C.push_back(e3);V_D.push_back(e4);}
  if (type == "F") print_TSIL_F(main_output, e1+e2+e3+e4+e5);
  
  main_output << "\n";
  
  }
  
  
  
  main_output << "}"  << endl;
  main_output << "\n";
  main_output << "\n";
  // init function
  
main_output << "  int init(Data data) \n"
<<"{\n"
<< "mw= data.M_w, mz = data.M_z ,ma = data.M_a, mc = data.M_chi, g2=data.g1;\n"
//<<"mw= 80.385L, mz = 91.1876L ,ma = 0.1, g2=0.65L;\n"
<<"tW = acos(mw/mz);\n"
<<"sw2 = TSIL_POW(sin(tW),2),S2TW=TSIL_POW(sin(tW),2), cw2 = TSIL_POW(cos(tW),2);\n"
<<"sw = sin(tW), cw = cos(tW);\n"
<<"C = TSIL_POW(g2,2)/(16.0L*PI*PI);\n"
<<"i=Power(-1,0.5);\n"
<<"Pi=PI;\n"
<<"return 0;\n"
<<"}\n";
  
  
  
  
  
  
  ifstream infile4("generator/functions.txt");
  
  string content2 = "";
  int ii;


  for(ii=0 ; infile4.eof()!=true ; ii++) // get content of infile
  {
  content2 += infile4.get();
  }
  ii--;
  content2.erase(content2.end()-1);     // erase last character

  cout << ii << " characters read...\n";
  infile4.close();

  main_output << content2;                 // output

  
  
  
  
  
  
  
  // MAIN FUNCTION CALL
  
  
  
main_output << "void Self_energy::init_tsil(Data data)\n"
<< "{\n"
<< "TSIL_REAL M=data.M_chi, qq= data.Q;\n"
<< "TSIL_REAL s=pow(M,2);\n"
<< "init(data);\n"
<< "DoTSIL(s,qq);\n"
<< "}\n"
<< "\n"
<< "\n"
<< "void Self_energy::run_tsil (Data &data) \n"
<< "{\n"
<< "TSIL_REAL M=data.M_chi;\n"
<< "p=M;\n"
<< "init_tsil(data);\n"
<< endl;

vector<std::string> particle_names_short = particle_names;


  sort(particle_names_short.begin(),particle_names_short.end());
  particle_names_short.erase( unique( particle_names_short.begin(), particle_names_short.end() ), particle_names_short.end() );




for (int i=0;i<particle_names_short.size();i++)
{
  main_output<< "TSIL_COMPLEXCPP SE_"<<particle_names_short[i]<<" = ";
  string particle_name_tmp = particle_names_short[i];
  for (int d = 0; d<nd;d++)
  {
  if (particle_names[d] == particle_name_tmp)
  {
  main_output<< " + diagram"<<"_"<< particle_name_tmp << "_" << tags[d] << "()";
  tag = tags[d];
  }
  
  
  }
  main_output << ";" << endl;
  main_output << "cout << \"Self energy of particle "<< particle_name_tmp << " = \" << SE_"<<particle_names_short[i]<<" << endl;"<<endl;
  
  main_output << "data.SE_"<< (i+1) << " = " << "real(SE_"<<particle_names_short[i]<<");"<<endl;
  
}

  main_output << "}" << endl;
  


  main_output.close();
  
  return 0;

  
  
}
