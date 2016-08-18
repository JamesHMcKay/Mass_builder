// to compile g++ -o stage_1 stage_1.cpp


/* goal: produce list of all possible permutations for objects F(x,y,z,u,v)
where (x,y,z,u,v) are in the set (a,b,c,d,e, ... ), somehow avoid repetitions
in this process */


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




// create a mathematica object that is the product of two basis integrals
void print_product(ofstream &myfile,string name_1,string name_2,string SEn="SEn")
{


  string elements1,elements2;
  stringstream _e1,_e2,_e3,_e4,_e5, _type1;
  string e1,e2,e3,e4,e5, type1;
  _type1 << name_1[0];
  _type1 >> type1;
  _e1 << name_1[1];_e1 >> e1;
  _e2 << name_1[2];_e2 >> e2;
  _e3 << name_1[3];_e3 >> e3;
  _e4 << name_1[4];_e4 >> e4;
  _e5 << name_1[5];_e5 >> e5;
  
  stringstream _f1,_f2,_f3,_f4,_f5, _type2;
  string f1,f2,f3,f4,f5, type2;
  _type2 << name_2[0];
  _type2 >> type2;
    
  _f1 << name_2[1];_f1 >> f1;
  _f2 << name_2[2];_f2 >> f2;
  _f3 << name_2[3];_f3 >> f3;
  _f4 << name_2[4];_f4 >> f4;
  _f5 << name_2[5];_f5 >> f5;


  if (type1=="A") elements1 = e1;
  if (type1=="B") elements1 = e1 + e2;
  if (type1=="J") elements1 = e1 + e2 + e3;
  if (type1=="T") elements1 = e1 + e2 + e3;
  if (type1=="K") elements1 = e1 + e2 + e3;
  if (type1=="V") elements1 = e1 + e2 + e3 + e4;

  if (type2=="A") elements2 = f1;
  if (type2=="B") elements2 = f1 + f2;
  if (type2=="J") elements2 = f1 + f2 + f3;
  if (type2=="T") elements2 = f1 + f2 + f3;
  if (type2=="K") elements2 = f1 + f2 + f3;
  if (type2=="V") elements2 = f1 + f2 + f3 + f4;



  myfile << type1 << elements1 << type2 << elements2 << " = ";
  if (type1=="A") myfile << "TAI[4, 0, {1, m" << elements1[0] << "}]";
  if (type1=="B") myfile << "TBI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements1[0] << "}, {1, m" << elements1[1] << "}}]";
  if (type1=="J") myfile << "TJI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements1[0] << "}, {1, m" << elements1[1] << "}, {1, m" << elements1[2] << "}}]";
  if (type1=="T") myfile << "TJI[4, Pair[Momentum[p],Momentum[p]], {{2, m" << elements1[0] << "}, {1, m" << elements1[1] << "}, {1, m" << elements1[2] << "}}]";
  if (type1=="K") myfile << "TJI[4, 0, {{1, m" << elements1[0] << "}, {1, m" << elements1[1] << "}, {1, m" << elements1[2] << "}}]";
  if (type1=="V") myfile << "TVI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements1[0] << "}, {1, m" << elements1[1] << "}, {1, m" << elements1[2] << "}, {1, m" << elements1[3] << "}}]";
  myfile << " * ";
  if (type2=="A") myfile << "TAI[4, 0, {1, m" << elements2[0] << "}];";
  if (type2=="B") myfile << "TBI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements2[0] << "}, {1, m" << elements2[1] << "}}];";
  if (type2=="J") myfile << "TJI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements2[0] << "}, {1, m" << elements2[1] << "}, {1, m" << elements2[2] << "}}];";
  if (type2=="T") myfile << "TJI[4, Pair[Momentum[p],Momentum[p]], {{2, m" << elements2[0] << "}, {1, m" << elements2[1] << "}, {1, m" << elements2[2] << "}}];";
  if (type2=="K") myfile << "TJI[4, 0, {{1, m" << elements2[0] << "}, {1, m" << elements2[1] << "}, {1, m" << elements2[2] << "}}];";
  if (type2=="V") myfile << "TVI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements2[0] << "}, {1, m" << elements2[1] << "}, {1, m" << elements2[2] << "}, {1, m" << elements2[3] << "}}];";
  myfile << endl;
  if ((type1 == type2) && (elements1==elements2)) myfile << "C"<< type1 << elements1 << type2 << elements2 << " = Coefficient["<<SEn<<","<< type1 << elements1 << ", 2];" << endl;
  else myfile << "C"<< type1 << elements1 << type2 << elements2 << " = - (1/2)* Coefficient["<<SEn<<","<< type1 << elements1 << type2 << elements2 << ", 1];" << endl;

}


void print_A(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "A"<< elements << " = " << "TAI[4, 0, {1, m" << elements[0] << "}];" << endl;
myfile << "CA"<< elements << " = Coefficient["<<SEn<<", A" << elements << ", 1];" << endl;
}

void print_B(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "B"<< elements << " = " << "TBI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements[0] << "}, {1, m" << elements[1] << "}}];" << endl;
myfile << "CB"<< elements << " = Coefficient["<<SEn<<", B" << elements << ", 1];" << endl;
}

void print_V(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "V"<< elements << " = " << "TVI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements[0] << "}, {1, m" << elements[1] << "}, {1, m" << elements[2] << "}, {1, m" << elements[3] << "}}];" << endl;
myfile << "CV"<< elements << " = Coefficient["<<SEn<<", V" << elements << ", 1];" << endl;
}

void print_T(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "T"<< elements << " = " << "TJI[4, Pair[Momentum[p],Momentum[p]], {{2, m" << elements[0] << "}, {1, m" << elements[1] << "}, {1, m" << elements[2] << "}}];" << endl;
myfile << "CT"<< elements << " = Coefficient["<<SEn<<", T" << elements << ", 1];" << endl;
}

void print_J(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "J"<< elements << " = " << "TJI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements[0] << "}, {1, m" << elements[1] << "}, {1, m" << elements[2] << "}}];" << endl;
myfile << "CJ"<< elements << " = Coefficient["<<SEn<<", J" << elements << ", 1];" << endl;
}

void print_K(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "K"<< elements << " = " << "TJI[4, 0, {{1, m" << elements[0] << "}, {1, m" << elements[1] << "}, {1, m" << elements[2] << "}}];" << endl;
myfile << "CK"<< elements << " = Coefficient["<<SEn<<", K" << elements << ", 1];" << endl;
}


void print_F(ofstream &myfile, string elements,string SEn = "SEn")
{
if (SEn == "SEn") myfile << "F"<< elements << " = " << "TFI[4, Pair[Momentum[p],Momentum[p]], {{1, m" << elements[0] << "}, {1, m" << elements[1] << "}, {1, m" << elements[2] << "}, {1, m" << elements[3] << "}, {1, m" << elements[4] << "}}];" << endl;
myfile << "CF"<< elements << " = Coefficient["<<SEn<<", F" << elements << ", 1];" << endl;
}



// TODO: we only need one list input here, and then to just duplicate the vector so we have it 5 times over so list input is just
// a list of the mass terms we require
void get_data(vector<std::string> &A, vector<std::string> &B , vector<std::string> &C , vector<std::string> &D, vector<std::string> &E, int &n)
{

std::ifstream input("list_input.txt");
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
A.resize(n);
B.resize(n);
C.resize(n);
D.resize(n);
E.resize(n);


n=0;
std::ifstream input2("list_input.txt");
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  iss2>> A[n] >> B[n] >> C[n] >> D[n] >> E[n];
    n=n+1;
 }




}



void get_data2(vector<std::string> &A,int &n)
{

std::ifstream input("names_updated.txt");
std::string line;
while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      n=n+1;
   }
  
A.resize(n);


n=0;
std::ifstream input2("names_updated.txt");
std::string line2;
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  iss2>> A[n];
    n=n+1;
 }




}






bool calc_diagram(string diagram,string particle)
{

  bool verbose=0;

  cout << "calculating diagram " << diagram << " for particle " << particle << endl;

  const char *ext = ".txt";



  
  time_t t = time(0);   // get time now to print into generated files
  struct tm * now = localtime( & t );
  
  
  
  
  string tag = particle + "_" + diagram;//"chi0_13";
  
  ofstream tag_out;
  tag_out.open ("tag.txt"); // store current tag = particle name + diagram number
  
  tag_out << tag << endl;
  
  tag_out.close();
  
  system("touch names_updated.txt");
  system("touch names_updated_temp.txt");
  
  int n = 0;
  vector<std::string> output0, output1,output2;
  vector<double> dup;
  vector<std::string> A;
  vector<std::string> B;
  vector<std::string> C;
  vector<std::string> D;
  vector<std::string> E;
  
  
  get_data(A, B , C ,D ,E, n);
  
  
  string s_cwd(getcwd(NULL,0));
  
  
  
  output0.resize(pow(n,5));
  
  
  ofstream myfile;
  myfile.open ("stage_3.m");


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
  <<"chi0 = InsertFields[t12, {F[6]} -> {F[6]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  <<"chi1 = InsertFields[t12, {F[5]} -> {F[5]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  <<"subdiags0 =   DiagramExtract["<<particle<<", "<<diagram<<"] (*57, 58, 59, 60, 61, 62];*)\n"
  <<"Export[\""<<s_cwd<<"/current_diagram.pdf\",Paint[subdiags0]];\n"  // print the FA diagram to pdf in local directory
  <<"amp0 := FCFAConvert[CreateFeynAmp[subdiags0], IncomingMomenta -> {p}, OutgoingMomenta -> {p}, LoopMomenta -> {k1, k2} ,UndoChiralSplittings -> True,(*TransversePolarizationVectors\[Rule] {k1},*)DropSumOver -> True, List -> False, ChangeDimension -> 4] // Contract\n"
  <<"amp0 = amp0 /. MajoranaSpinor[p, mc] -> 1 /.Spinor[Momentum[p], mc, 1] -> 1;\n"
  <<"SetOptions[Eps, Dimension -> D];\n"
  <<"fullamp0 = (amp0) // DiracSimplify // FCMultiLoopTID[#, {k1, k2}] & //DiracSimplify;\n"
  <<"tfiamp0 = fullamp0 // ToTFI[#, k1, k2, p] & // ChangeDimension[#, 4] &;\n"
  <<"Print[tfiamp0]\n"
  <<"SEn = FullSimplify[TarcerRecurse[tfiamp0] /. D -> 4 /.MajoranaSpinor[p, mc] -> 1] /. Spinor[Momentum[p], mc, 1] -> 1;\n"
  <<"DumpSave[\"stage_3.mx\", SEn];\n"
  <<"Print[SEn]"<< endl;


  ofstream myfile_names;
  myfile_names.open ("names.txt");
  
  vector<std::string> Bases;
  int nb = 0;
  
  int k=0;
  
  // want to do for loop for each element, but we have arbitrary number of elements
  
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  print_A(myfile,A[i1]);
  myfile_names << "A"<< A[i1] << endl;
  Bases.push_back("A"+A[i1]);nb = nb+1;
  for (int i2 = 0; i2 < n ; i2++)
  {
  for (int i3 = 0; i3 < n ; i3++)
  {
  for (int i4 = 0; i4 < n ; i4++)
  {
  print_V(myfile , A[i1]+B[i2]+C[i3]+D[i4]);
  myfile_names << "V"<< A[i1]+B[i2]+C[i3]+D[i4] << endl;
  Bases.push_back("V"+A[i1]+B[i2]+C[i3]+D[i4]);nb = nb+1;
  
  for (int i5 = 0; i5 < n ; i5++)
  {
  print_F(myfile, A[i1]+B[i2]+C[i3]+D[i4]+E[i5] );
  myfile_names << "F"<< A[i1]+B[i2]+C[i3]+D[i4]+E[i5] << endl;
  Bases.push_back("F"+A[i1]+B[i2]+C[i3]+D[i4]+E[i5]);nb = nb+1;
  
  k=k+1;
  }}}}}
  
  
  
  // deal with symmetries appropriately
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  print_B(myfile, A[i1]+A[i2]);
  myfile_names << "B"<< A[i1]+B[i2] << endl;
  Bases.push_back("B"+A[i1]+B[i2]);nb = nb+1;
  
  }}
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  print_J(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "J"<< A[i1]+B[i2]+C[i3] << endl;
  Bases.push_back("J"+A[i1]+B[i2]+C[i3]);nb = nb+1;
  print_K(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "K"<< A[i1]+B[i2]+C[i3] << endl;
  Bases.push_back("K"+A[i1]+B[i2]+C[i3]);nb = nb+1;
  
  }}}
  
    
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = 0; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  
  print_T(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "T"<< A[i1]+B[i2]+C[i3] << endl;
  Bases.push_back("T"+A[i1]+B[i2]+C[i3]);nb = nb+1;
  
  
  }
  }
  }
  
  
  
  myfile << "Export[\""<<s_cwd<<"/output.txt\", {" << endl;
  
  
  
  for (int i=0;i<nb-1;i++)
  {
    myfile << "{\" TSIL_COMPLEXCPP C"<<Bases[i]<<" =\", CForm[C"<<Bases[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  }
  
  myfile << "{\" TSIL_COMPLEXCPP C"<<Bases[nb-1]<<" =\", CForm[C"<<Bases[nb-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  myfile << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  myfile << "Print[\"Completed\"]"<<endl;
  
  
  
  myfile.close();
  
  myfile_names.close();
  

  #ifdef RUN_ALL
  system("chmod +x stage_3.m ");
  if (verbose) system("./stage_3.m");
  else system("./stage_3.m  >/dev/null");
  system("./stage_4.sh");
  #endif
  
  
  
  
  
  //////////////////////////////////////////////////////////////////
  //  Now determine create substatially more compact Mathematica routine
  //  using only the required basis integrals and non-zero coefficients
  //  then determine what cross-terms we still need to account for
  
  
  
  n = 0;
  vector<std::string> integrals;
  
  
  get_data2(integrals, n);
  
  
  string prev = "\"stage_3.mx\"";
  ofstream myfile_stage6;
  myfile_stage6.open ("stage_6.m");
  
  myfile_stage6 << "#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script\n"
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
  <<"Get[" << prev << "]\n"
  <<"Print[SEn]"<< endl;
  
  
  
  
  
  // integrals is a vector of length n containing all the required integrals
  // need to extract each element and catagorise them
  
  for (int i = 0; i<n;i++)
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
  
  if (type == "A") print_A(myfile_stage6,e1);
  if (type == "B") print_B(myfile_stage6, e1+e2);
  if (type == "J") print_J(myfile_stage6, e1+e2+e3);
  if (type == "T") print_T(myfile_stage6, e1+e2+e3);
  if (type == "K") print_K(myfile_stage6, e1+e2+e3);
  if (type == "V") print_V(myfile_stage6, e1+e2+e3+e4);
  if (type == "F") print_F(myfile_stage6, e1+e2+e3+e4+e5);
  }

  
  myfile_stage6 << "SEnTrial = ";
  
  for (int i = 0; i<n;i++)
  {
  
  myfile_stage6 << " + "<<integrals[i]<< " * C"<<integrals[i];
    
  }
  myfile_stage6 << ";"<<endl;
  
  myfile_stage6 << "diff = FullSimplify[SEn-SEnTrial]"<<endl;
  myfile_stage6 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage6 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage6 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  myfile_stage6 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage6 << "Print[\" --------------------------------------- \"]" <<endl;
  
  myfile_stage6 << "Print[diff]"<<endl;
  
  
  
  

  
  

    
  for (int i = 0; i<n;i++)
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
  
  if (type == "A") print_A(myfile_stage6,e1,"diff");
  if (type == "B") print_B(myfile_stage6, e1+e2,"diff");
  if (type == "J") print_J(myfile_stage6, e1+e2+e3,"diff");
  if (type == "T") print_T(myfile_stage6, e1+e2+e3,"diff");
  if (type == "K") print_K(myfile_stage6, e1+e2+e3,"diff");
  if (type == "V") print_V(myfile_stage6, e1+e2+e3+e4,"diff");
  if (type == "F") print_F(myfile_stage6, e1+e2+e3+e4+e5,"diff");
  }
  
  
  myfile_stage6 << "Export[\""<<s_cwd<<"/output.txt\", {" << endl;
  
    for (int i=0;i<n-1;i++)
  {
    myfile_stage6 << "{\" TSIL_COMPLEXCPP C"<<integrals[i]<<" =\", CForm[C"<<integrals[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  }
  
  myfile_stage6 << "{\" TSIL_COMPLEXCPP C"<<integrals[n-1]<<" =\", CForm[C"<<integrals[n-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  myfile_stage6 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  myfile_stage6 << "Print[\"Completed\"]"<<endl;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  myfile_stage6.close();
  
  
  #ifdef RUN_ALL
  system("chmod +x stage_6.m ");
  if(verbose) system("./stage_6.m ");
  else system("./stage_6.m  >/dev/null ");
  system("chmod u+x stage_8.sh");
  system("./stage_8.sh ");
  #endif
  
  
  //////////////////////////////////////////////////////////////////
  // of the terms that are remaining we now deal with possible combinations, the updated file for objects of interest
  // is names_updated.txt, this should contain all integrals for which we need to consider combinations
  
  int m = n;
  n = 0;
  vector<std::string> products;
  get_data2(products, n);
  
  // need to consider all permuations of each element of combos, for which rank < 5
  // ranks are scored as ->  J,T,K = 3, A = 1, B = 2, V = 4, F = 5 (should not come up here)
  // produce new objects, for example the object BacAc = TAI [ ... ] * TBI [ ... ] and the corresponding coefficient check
  // and print out statement in CForm
  
  
  
  ofstream myfile_stage8;
  myfile_stage8.open ("stage_8.m");
  
  
  
  prev = "\"stage_3.mx\"";
  myfile_stage8 << "#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script\n"
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
  <<"Get[" << prev << "]\n"
  <<"Print[SEn]"<< endl;

    
  for (int i = 0; i<m;i++)
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
  
  if (type == "A") print_A(myfile_stage8,e1);
  if (type == "B") print_B(myfile_stage8, e1+e2);
  if (type == "J") print_J(myfile_stage8, e1+e2+e3);
  if (type == "T") print_T(myfile_stage8, e1+e2+e3);
  if (type == "K") print_K(myfile_stage8, e1+e2+e3);
  if (type == "V") print_V(myfile_stage8, e1+e2+e3+e4);
  if (type == "F") print_F(myfile_stage8, e1+e2+e3+e4+e5);
  }
  
  
  
  
  for (int i = 0; i<n;i++)
  {

  
  for (int j = 0; j<n;j++)
  {
  
  
  print_product(myfile_stage8,products[i],products[j]);
  }
  }
  myfile_stage8 << "SEnTrial = ";
  for (int i = 0; i<m;i++)
  {
  myfile_stage8 << " + "<<integrals[i]<< " * C"<<integrals[i];
  }
  
  
  
  
  for (int i = 0; i<n;i++)
  {
  for (int j = 0; j<n;j++)
  {
  myfile_stage8 << " + "<< products[i] << products[j] << " * C"<<products[i]<< products[j];
  }
  }
  myfile_stage8 << ";"<<endl;

  
  myfile_stage8 << "diff = FullSimplify[SEn-SEnTrial]"<<endl;
  myfile_stage8 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage8 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage8 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  myfile_stage8 << "Print[\" --------------------------------------- \"]" <<endl;
  myfile_stage8 << "Print[\" --------------------------------------- \"]" <<endl;
  
  myfile_stage8 << "Print[diff]"<<endl;
  myfile_stage8 << "Export[\"result.txt\",diff]"<<endl;
  
  // print out coefficients of products
  
  myfile_stage8 << "Export[\""<<s_cwd<<"/output_products.txt\", {" << endl;
  
  for (int i = 0; i<n;i++)
  {
  for (int j = 0; j<n;j++)
  {
    
    if ((i==n-1) && (j==n-1)){ myfile_stage8 << "{\" TSIL_COMPLEXCPP C"<<products[i] << products[j] <<" =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;}
    else {myfile_stage8 << "{\" TSIL_COMPLEXCPP C"<<products[i] << products[j] <<" =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
    
  }
  }
  
  myfile_stage8 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  // print out in a form useful for determing products later on
  
  myfile_stage8 << "Export[\""<<s_cwd<<"/output_products_2.txt\", {" << endl;
  
  for (int i = 0; i<n;i++)
  {
  for (int j = 0; j<n;j++)
  {
    
    if ((i==n-1) && (j==n-1)){ myfile_stage8 << "{\" "<<products[i] << "*" << products[j] << "*C"<<products[i] << products[j] <<"        =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
    else {myfile_stage8 << "{\" "<<products[i] << "*" << products[j] << "*C"<<products[i] << products[j] <<"           =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
    
  }
  }
  
  myfile_stage8 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;


  myfile_stage8.close();
  
  #ifdef RUN_ALL
  system("chmod +x stage_8.m ");
  if(verbose) system("./stage_8.m");
  else system("./stage_8.m  >/dev/null ");
  system("chmod u+x stage_9.sh");
  system("./stage_9.sh");
  #endif
  
  std::ifstream file("result.txt");
  std::string str;
  std::string result;
  std::getline(file, str);
  result += str;
  
 
  // need to check if the result of diff is zero, if not then throw an error
  
  bool success = 0;
  if (result == "0") {cout << "Successful!!!" << endl, success = 1;}
  else { cout << "Something has gone terribly wrong, check the symmetries are being accounted for properly when writting the basis integrals \n" <<\
  "(and thus no double up of terms) and that all cross terms are being considered." << endl; success = 0;}
  
  
  
  
  
  // need a list of basis integrals that are required and a list of the coefficients
  // need the string SEn = " ... " for use in writting the C++ code
  
  
  
  const char* basis_integrals_tmp = "output/basis_integrals_"; // vector containing file names
  string c_basis_integrals = basis_integrals_tmp + tag + ext;
  const char *basis_integrals = c_basis_integrals.c_str();


  ofstream output_1;
  output_1.open (basis_integrals);

  for (int i = 0; i<m;i++)
  {
  output_1 << integrals[i] << endl;
  }

  output_1.close();
  
  ofstream diagram_list;
  diagram_list.open( "output/diagrams.txt", ios::out | ios::app );
  diagram_list << particle << " " << diagram << endl;
  diagram_list.close();
  
  
  
  // write full list of basis integral * coefficient for non-zero terms
  
  
  
  const char* summation_tmp = "output/summation_"; // vector containing file names
  string c_summation = summation_tmp + tag + ext;
  const char *summation = c_summation.c_str();

  ofstream summation_out;
  
  summation_out.open (summation);

  summation_out << "return ";
  for (int i = 0; i<m;i++)
  {
  summation_out  << " + "<<integrals[i]<< " * C"<<integrals[i];
  }
  
  
  int np = 0;
  vector<std::string> products_reduced;
  get_data2(products_reduced, np);
  
  
  for (int i = 0; i<np;i++)
  {
  
  
  summation_out << " + "<< products_reduced[i];
  
  }

  
  
  
  
  summation_out << ";"<<endl;
  summation_out.close();
  
  
  
  
  
  
  
  return success;
}


// main routine to manage a diagram by diagram procedure


int main(int argc, char *argv[])
{

string diagram = "";
string particle = "";
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
  std::string particles [100];
  std::string diagrams [100]; int i=0;
  std::ifstream input(argv[2]);
  std::string line;
  while(getline(input, line))
  {
        if (!line.length() || line[0] == '#')
           continue;
        std::istringstream iss(line);
        iss>> particles[i] >> diagrams[i];
        i=i+1;
  }
// run over all entries in the input file
  for (int k=0;k<i;k++)
  {
  
  calc_diagram(diagrams[k],particles[k]);
  }

}
else
{
diagram = argv[2];
particle = argv[1];
calc_diagram(diagram,particle);
}


}


}
