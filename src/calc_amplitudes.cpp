/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 28/09/16
James McKay

--- calc_amplitudes.cpp ---

generate mathematica scripts to compute all coefficients and find required basis integrals
also includes functions for sending FeynArts diagrams directory to pdf in the folder FAdiagrams

run ./mass_builder -g <model_name> after to generate code and ./mass_builder -e input.txt to evaluate

*/

#include "calc_amplitudes.hpp"

#define RUN_ALL

using namespace std;
using namespace utils;

// create a mathematica object that is the product of two basis integrals

bool verbose=1;
int loop_level = 2;




bool Calc_amplitudes::calc_diagram(string diagram,string particle,string model)
{
  bool success=0;
  bool sum_integrals=1;
  

  cout << "calculating diagram " << diagram << " for particle " << particle << " in model " << model << endl;


  const char *ext = ".txt";
  string underscore = "_";
  string blank = "";

  string particle_full = particle;
 
  particle =  part_name_simple(particle_full);
  // edit particle name to a safe string
  
  
  
  
  string tag = particle + "_" + diagram;//"chi0_13";
  
  ofstream tag_out;
  tag_out.open ("output/tag.txt"); // store current tag = particle name + diagram number
  
  tag_out << tag << endl;
  
  tag_out.close();
  
  ofstream model_out;
  model_out.open ("output/model.txt"); // store current tag = particle name + diagram number
  
  model_out << model << endl;
  
  model_out.close();
  
  
  
  system("touch output/names_updated.txt");
  system("touch output/names_updated_temp.txt");
  
  int n = 0;
  vector<std::string> output0, output1,output2;
  vector<double> dup;
  vector<std::string> A;
  
  const char* file_masses_tmp = "models/";
  string c_file_masses = file_masses_tmp + model + "/masses" + ext;
  const char *file_masses = c_file_masses.c_str();
  
  
  get_data(A, n,file_masses);
  
  
  string s_cwd(getcwd(NULL,0));
  
  
  
  output0.resize(pow(n,5));
  
  
  ofstream myfile;
  myfile.open ("output/stage_3.m");



  utils::print_math_header(myfile);
  utils::print_math_body(myfile,loop_level,particle_full,diagram, model,s_cwd);
 
  myfile<<"Print[tfiamp0]\n"
  <<"SEn = FullSimplify[TarcerRecurse[tfiamp0] /. D -> 4 /.MajoranaSpinor[p, mc] -> 1] /. Spinor[Momentum[p], mc, 1] -> 1;\n"
  <<"DumpSave[\""<<s_cwd<<"/output/stage_3.mx\", SEn];\n"
  <<"Print[\"----------- The self energy is ---------- = \"]\n"
  <<"Print[SEn]\n"
  <<"Print[\"-------------------- = \"]"<< endl;


  
  
  
  


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
  
  plan to make handling of basis integrals more abstract to deal with complicated mass names
  
  we create a struct for each basis integral that is given a simplifed name, this struct contains the mass
  in each location of the integral, and the type of basis integral
  
  simplifed name: FMaMaMaMzMw is just stupid, it should be Faaazw, so we need to assign a unique single character
  name to each mass, this can be set at runtime by the user in the second column of "masses", if it is not set
  it will be set automatically from the alphabet.
  
  alternatively, we have a class for each type of basis integral which is initialised with the required masses
  the class then determines all possible non degenerate basis integrals and creates these in vectors (or it could
  make a set of structs like that described above)
  
  
  tricky part is identifying the non-zero elements after a run through mathematica, what we get back is
  at best a list of the non-zero class names, for example Faaazw.  It would be nicer if somehow the struct
  for each basis integral held a string called "coeff", which was automatically updated after a run through
  mathematica.  Then all those with zeros are deallocated.
  
  idea to solve above problem: open mathematica output, and directly read in the following
  "Faaazw COEFFICIENT"
  
  so someting like
  
  open(math_out.txt)
  
  line by line:  string class_name << Faaazw
  string COEFF << COEFFICIENT
  
  
  util_func( class_name, COEFF) { writes COEFF into class given by class_name} // tricky -> use maps
  so for each mass combination we create a class with a 
  
  
  
  
  
  */
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  ofstream myfile_names;
  myfile_names.open ("output/names.txt");
  
  vector<std::string> Bases;
  int nb = 0;
  
  int k=0;
  
  
  
  
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
  print_V(myfile , A[i1]+A[i2]+A[i3]+A[i4]);
  myfile_names << "V"<< A[i1]+A[i2]+A[i3]+A[i4] << endl;
  Bases.push_back("V"+A[i1]+A[i2]+A[i3]+A[i4]);nb = nb+1;
  
  for (int i5 = 0; i5 < n ; i5++)
  {
  print_F(myfile, A[i1]+A[i2]+A[i3]+A[i4]+A[i5] );
  myfile_names << "F"<< A[i1]+A[i2]+A[i3]+A[i4]+A[i5] << endl;
  Bases.push_back("F"+A[i1]+A[i2]+A[i3]+A[i4]+A[i5]);nb = nb+1;
  
  k=k+1;
  }}}}}
  
  
  
  // deal with symmetries appropriately
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  print_B(myfile, A[i1]+A[i2]);
  myfile_names << "B"<< A[i1]+A[i2] << endl;
  Bases.push_back("B"+A[i1]+A[i2]);nb = nb+1;
  
  }}
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  print_J(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "J"<< A[i1]+A[i2]+A[i3] << endl;
  Bases.push_back("J"+A[i1]+A[i2]+A[i3]);nb = nb+1;
  print_K(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "K"<< A[i1]+A[i2]+A[i3] << endl;
  Bases.push_back("K"+A[i1]+A[i2]+A[i3]);nb = nb+1;
  
  }}}
  
    
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = 0; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  
  print_T(myfile, A[i1]+A[i2]+A[i3]);
  myfile_names << "T"<< A[i1]+A[i2]+A[i3] << endl;
  Bases.push_back("T"+A[i1]+A[i2]+A[i3]);nb = nb+1;
  
  
  }
  }
  }
  
  myfile_names.close();
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  myfile << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;
  
  
  
  for (int i=0;i<nb-1;i++)
  {
    myfile << "{\" TSIL_COMPLEXCPP C"<<Bases[i]<<" =\", CForm[C"<<Bases[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  }
  
  myfile << "{\" TSIL_COMPLEXCPP C"<<Bases[nb-1]<<" =\", CForm[C"<<Bases[nb-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;
  myfile << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  myfile << "Print[\"Completed\"]"<<endl;
  
  
  
  myfile.close();
  
  
  

  #ifdef RUN_ALL
  system("chmod +x output/stage_3.m ");
  if (verbose) system("./output/stage_3.m");
  else system("./output/stage_3.m  >/dev/null");
  system("chmod u+x scripts/script_1.sh");
  system("./scripts/script_1.sh");
  #endif
  
  
  
  
  
  //////////////////////////////////////////////////////////////////
  //  Now create substatially more compact Mathematica routine
  //  using only the required basis integrals and non-zero coefficients
  //  then determine what cross-terms we still need to account for
  
  
  
  n = 0;
  vector<std::string> integrals;
  
  
  const char* file_integrals_tmp = "output/names_updated";
  string c_file_integrals = file_integrals_tmp + blank + ext;
  const char *file_integrals = c_file_integrals.c_str();
  
  
  get_data(integrals, n,file_integrals);
  if (n==0) {integrals = Bases; n = Bases.size(); sum_integrals = 0;}
  
  
  string prev = "stage_3.mx\"";
  ofstream myfile_stage6;
  myfile_stage6.open ("output/stage_6.m");
  
  utils::print_math_header(myfile_stage6);
  myfile_stage6<<"Get[\"" << s_cwd <<"/output/"<< prev << "]\n"
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

  
  myfile_stage6 << "SEnTrial = 0 ";
  
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
  
  myfile_stage6 << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;

  
  

    
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
  
  
  myfile_stage6 << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;
  if (n>0)
  {
    for (int i=0;i<n-1;i++)
  {
    myfile_stage6 << "{\" TSIL_COMPLEXCPP C"<<integrals[i]<<" =\", CForm[C"<<integrals[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
  }
  
  myfile_stage6 << "{\" TSIL_COMPLEXCPP C"<<integrals[n-1]<<" =\", CForm[C"<<integrals[n-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;
  myfile_stage6 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  }
  
  myfile_stage6 << "Print[\"Completed\"]"<<endl;
  
  
  
  
  
  
  
  
  
  
  myfile_stage6.close();
  
  
  #ifdef RUN_ALL
  system("chmod +x output/stage_6.m ");
  if(verbose) system("./output/stage_6.m ");
  else system("./output/stage_6.m  >/dev/null ");
  system("chmod u+x scripts/script_2.sh");
  system("./scripts/script_2.sh ");
  #endif
  
  int m = n;
  n = 0;
  
//  if (check_done()) {success = 1;}
 // else
//  {
  
  
  //////////////////////////////////////////////////////////////////
  // of the terms that are remaining we now deal with possible combinations, the updated file for objects of interest
  // is names_updated.txt, this should contain all integrals for which we need to consider combinations
  

  vector<std::string> products;
  get_data(products, n,file_integrals);
  
  if (n == 0){
  // just use Bases instead but remove all F terms!
  
  for (unsigned int i=0;i<Bases.size();i++)
  {
  
  string nameB = Bases[i];
  stringstream _Btype;
  string Btype;
  _Btype << nameB[0];
  _Btype >> Btype;
  
  if (Btype != "F") { n=n+1; products.push_back(Bases[i]);}
  
  
  }
  }
  
  
  // need to consider all permuations of each element of combos, for which rank < 5
  // ranks are scored as ->  J,T,K = 3, A = 1, B = 2, V = 4, F = 5 (should not come up here)
  // produce new objects, for example the object BacAc = TAI [ ... ] * TBI [ ... ] and the corresponding coefficient check
  // and print out statement in CForm
  
  
  
  ofstream myfile_stage8;
  myfile_stage8.open ("output/stage_8.m");
  
  
  
  prev = "stage_3.mx\"";
  utils::print_math_header(myfile_stage8);
  myfile_stage8<<"Get[\"" << s_cwd <<"/output/"<< prev << "]\n"
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
  
  if (type == "A") print_A(myfile_stage8, e1);
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
  myfile_stage8 << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;
  
  // print out coefficients of products
  
  myfile_stage8 << "Export[\""<<s_cwd<<"/output/output_products.txt\", {" << endl;
  
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
  
  myfile_stage8 << "Export[\""<<s_cwd<<"/output/output_products_2.txt\", {" << endl;
  
  for (int i = 0; i<n;i++)
  {
  for (int j = 0; j<n;j++)
  {
    
    if ((i==n-1) && (j==n-1)){ myfile_stage8 << "{\" "<<products[i] << "*" << products[j] << "*C"<<products[i] << products[j] <<"        =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;}
    else {myfile_stage8 << "{\" "<<products[i] << "*" << products[j] << "*C"<<products[i] << products[j] <<"           =\", CForm[C"<<products[i] << products[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
    
  }
  }
  
  myfile_stage8 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;


  myfile_stage8.close();
  
  #ifdef RUN_ALL
  system("chmod +x output/stage_8.m ");
  if(verbose) system("./output/stage_8.m");
  else system("./output/stage_8.m  >/dev/null ");
  system("chmod u+x scripts/script_3.sh");
  system("./scripts/script_3.sh");
  #endif
 
  success = check_done();
  
 
  
  
  
  
  // need a list of basis integrals that are required and a list of the coefficients
  // need the string SEn = " ... " for use in writting the C++ code
  
  
  
  const char* basis_integrals_tmp = "/output/basis_integrals_"; // vector containing file names
  string c_basis_integrals = "models/" + model + basis_integrals_tmp + tag + ext;
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
  
  
  
  const char* summation_tmp = "/output/summation_"; // vector containing file names
  string c_summation = "models/" + model +summation_tmp + tag + ext;
  const char *summation = c_summation.c_str();

  ofstream summation_out;
  
  summation_out.open (summation);

  summation_out << "return ";
  for (int i = 0; i<m;i++)
  {
  if (sum_integrals != 0 ) summation_out  << " + "<<integrals[i]<< " * C"<<integrals[i];
  }
  
  
  int np = 0;
  vector<std::string> products_reduced;
  get_data(products_reduced, np,file_integrals);
  
  
  for (int i = 0; i<np;i++)
  {
  
  
  summation_out << " + "<< products_reduced[i];
  
  }

  
  
  
  
  summation_out << ";"<<endl;
  summation_out.close();
  
  
  
  
  
  
  
  return success;
}


void draw_all_diagrams(std::string particle, string model)
{

  



  
  
  string s_cwd(getcwd(NULL,0));
  
  
  
  ofstream myfile;
  myfile.open ("output/make_figures.m");


  utils::print_math_header(myfile);
  myfile<<"t12 = CreateTopologies["<<loop_level<<", 1 -> 1, ExcludeTopologies -> Internal];\n"
  <<"alldiags = InsertFields[t12, {"<<particle<<"} -> {"<<particle<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<s_cwd<<"/models/"<<model<<"/"<<model<<"\"];\n"
  <<"Export[\""<<s_cwd<<"/FA_diagrams/all_diagrams_"<<particle<<".pdf\",Paint[alldiags]];\n"  // print the FA diagram to pdf in local directory
  <<endl;

  #ifdef RUN_ALL
  system("chmod +x output/make_figures.m ");
  if (verbose) system("./output/make_figures.m");
  else system("./output/make_figures.m  >/dev/null");
  #endif
  
  
  
}




void draw_diagrams(vector<std::string> particles, vector<std::string> diagrams, int nd,string model)
{

  

  
  string s_cwd(getcwd(NULL,0));
  
  
  
  
  
  ofstream myfile;
  myfile.open ("output/make_figures.m");


 utils::print_math_header(myfile);
 
 myfile <<"t12 = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> Internal];\n"
  //<<"chi0 = InsertFields[t12, {F[6]} -> {F[6]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  //<<"chi1 = InsertFields[t12, {F[5]} -> {F[5]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \"MDM_tripletEWSB\"];\n"
  <<endl;
  
  
  vector<std::string> particle_names_short = particles;


  sort(particle_names_short.begin(),particle_names_short.end());
  particle_names_short.erase( unique( particle_names_short.begin(), particle_names_short.end() ), particle_names_short.end() );
  
  for (unsigned int i=0;i<particle_names_short.size();i++)
  {
  string particle_name_tmp = particle_names_short[i];
  myfile<<particle_name_tmp<<" = InsertFields[t12, {"<<particle_name_tmp<<"} -> {"<<particle_name_tmp<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<s_cwd<<"/models/"<<model<<"/"<<model<<"\"];\n";
  myfile<< "subdiags" << particle_name_tmp <<" =   DiagramExtract["<<particle_name_tmp;
  for (int d = 0; d<nd;d++)
  {
  if (particles[d] == particle_name_tmp)
  {
  
  myfile <<", "<<diagrams[d];
  //first =0;
  }
  }
  myfile <<"]\n";
  //first = 1;
  myfile <<"Export[\""<<s_cwd<<"/FA_diagrams/subset_diagrams_"<<particle_name_tmp<<".pdf\",Paint[subdiags"<<particle_name_tmp<<"]];\n"  // print the FA diagram to pdf in local directory
  <<endl;
  }

  #ifdef RUN_ALL
  system("chmod +x output/make_figures.m ");
  if (verbose) system("./output/make_figures.m");
  else system("./output/make_figures.m  >/dev/null");
  #endif
  
  
  
}


// main routine to manage a diagram by diagram procedure


void Calc_amplitudes::generate_figures(int argc, char *argv[])
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
if (argc==4) {cout <<"drawing all diagrams for particle "<<argv[2]<<endl;;draw_all_diagrams(argv[2],argv[3]);}
else {cout << "please specify a particle after -a to use this option"<< endl;}


}


}


}
