/*
 Mass Builder 
 
 James McKay
 Aug - Sep 2016
 Jan 2017
 
 --- calc_amplitudes.cpp ---
 
 This file contains the main routine for calling FeynArts, FeynCalc and TARCER via
 generated Mathematica scripts, followed by sorting this output and forming further output
 in the form of coefficient * basis_integral.
 
  This file also includes functions for sending FeynArts diagrams directory to pdf in the folder FAdiagrams
 */

#include "calc_amplitudes.hpp"

#define RUN_ALL
//#define DEBUG
using namespace std;
using namespace utils;
string s_cwd(getcwd(NULL,0));

void Calc_amplitudes::compute_amp(string prevb,string dimension)
{
  //  read in data from first run and reduce the basis of non-product terms //
  
  vector<std::string> output_string, coeff_new;
  int temp_int = 0; // this should always equal nb in this case
  
  const char* file_integrals2_tmp = "output/output";
  string c_file_integrals2 = file_integrals2_tmp + blank + ext;
  const char *file_integrals2 = c_file_integrals2.c_str();
  
  get_data(output_string, coeff_new, temp_int,file_integrals2, true);
  for (int i=0; i<temp_int; i++)
  {
    full_basis[output_string[i]].coefficient = coeff_new[i]; // set coefficients from first math run
  }
  
  reduced_basis  = remove_zeros(full_basis, full_basis_id); // remove integrals with zero coefficients
  reduced_basis_id = extract_keys(reduced_basis);
  
  nbr = reduced_basis_id.size();
  
  //  done updating basis map //
  
  // second run is an attempt to construct the trial amplitude using the reduced list of required integrals
  
  ofstream math_2;
  math_2.open ("output/math_2.m");
  
  utils::print_math_header(math_2);
  
  math_2<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n";
  
  print_math_basis(reduced_basis,math_2, "SEn",dimension);
  
  math_2 << "SEnTrial = 0 ";
  for (int i = 0; i<nbr;i++)
  {
    math_2 << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i]; // use the terms we no are non-zero
  }
  math_2 << ";"<<endl;
  
  math_2 << "difference = Simplify[SEn-SEnTrial]"<<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  math_2 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  math_2 << "Print[difference]"<<endl;
  math_2 << "Export[\""<<s_cwd<<"/output/result.txt\", difference]" << endl;
  
  print_math_basis(full_basis,math_2, "difference",dimension);

  
  math_2 << "Export[\""<<s_cwd<<"/output/output2.txt\", {" << endl;
  for (int i = 0; i < nb-1;i++)
  {
    math_2 << "{\""<<full_basis_id[i]<<" \", CForm[C"<<full_basis_id[i]<<" + C" <<full_basis_id[i] <<"2 /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  math_2 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C"<<full_basis_id[nb-1]<<" + C"<< full_basis_id[nb-1] <<"2 /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  math_2 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  math_2.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_2.m ");
  if(options.verbose) system("./output/math_2.m ");
  else system("./output/math_2.m  >/dev/null ");
#endif
  
  
  // part 2 complete //
  
  
  temp_int = 0;  // this should be equal to nb always
  const char* file_integrals4_tmp = "output/output2";
  string c_file_integrals4 = file_integrals4_tmp + blank + ext;
  const char *file_integrals4 = c_file_integrals4.c_str();
  vector<std::string> coeff_new_prod;
  prod_basis.clear();
  prod_id.clear();
  get_data(prod_id, coeff_new_prod, temp_int,file_integrals4,true);
  
  prod_basis = remove_type_F(full_basis,full_basis_id);
  
  if (temp_int == 0)  // this case should be redundant now
  {
    cout << "ERROR THIS CASE SHOULD BE REDUNDANT" << endl;
  }
  else
  {
    for (int i=0; i<temp_int; i++)
    {
      prod_basis[prod_id[i]].coefficient = coeff_new_prod[i];
    }
    prod_basis = remove_zeros(prod_basis, prod_id);
    
    prod_id = extract_keys(prod_basis);
    prod_basis = remove_type_F(prod_basis, prod_id);
    
    prod_id = extract_keys(prod_basis);
  }
  
  if (prod_basis.size()==0 && !check_done_quiet())
  {
    prod_basis = reduced_basis;//full_basis;
    prod_id = extract_keys(prod_basis);
    prod_basis = remove_type_F(prod_basis, prod_id);
    prod_id = extract_keys(prod_basis);
  }
  
  np = prod_id.size();
}

void Calc_amplitudes::make_full_trial(string prevb,string dimension,bool cform)
{

  ofstream math_3;
  math_3.open ("output/math_3.m");
  
  utils::print_math_header(math_3);
  math_3<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n";
  
  print_math_basis(reduced_basis,math_3,"SEn",dimension);
  print_math_basis(prod_basis,math_3,"SEn",dimension);
  print_math_products(prod_basis,math_3,"SEn",dimension);
  
  
  math_3 << "SEnTrial = 0 ";
  for (int i = 0; i<nbr;i++)
  {
    math_3 << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i];
  }
  
  for (int i = 0; i<np ; i++)
  {
    for (int j = 0; j<np ; j++)
    {
      math_3 << " + "<< prod_id[i] << prod_id[j] << " * C"<<prod_id[i]<< prod_id[j];
    }
  }
  math_3 << ";"<<endl;
  
  
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[\" The trial SE is:\"]" <<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[SEnTrial]"<<endl;
  
  math_3 << "diff = Simplify[SEn-SEnTrial]"<<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[diff]"<<endl;
  
  math_3 << "Export[\""<<s_cwd<<"/output/result.txt\", CForm[diff/. DiracGamma[Momentum[p]] -> p] ]" << endl;
  
  math_3<<"remainder = diff;"<<endl;
  
  math_3<<"DumpSave[\""<<s_cwd<<"/output/remainder.mx\", remainder];"<<endl;
  
  // print out coefficients of products -- this is for the purposes of passing as final output to TSIL
  
  math_3 << "Export[\""<<s_cwd<<"/output/output_products.txt\", {" << endl;
  
  
  products_map.clear();
  
  
  
  if (cform)
  {
    for (int i = 0; i<np;i++)
    {
      for (int j = 0; j<np;j++)
      {
        if ((i==np-1) && (j==np-1)){ math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C"<<prod_id[i] << prod_id[j] <<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;}
        else {math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C"<<prod_id[i] << prod_id[j] <<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;}
        
        // need to store information regarding the integrals involved in each product
        Bases_product product(prod_basis[prod_id[i]],prod_basis[prod_id[j]],prod_id[i],prod_id[j]);
        products_map[prod_id[i] + prod_id[j]] = product;
      
      }
    }
  }
  else
  {
    for (int i = 0; i<np;i++)
    {
      for (int j = 0; j<np;j++)
      {
        if ((i==np-1) && (j==np-1)){ math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",C"<<prod_id[i] << prod_id[j] <<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p, \"\"}" << endl;}
        else {math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",C"<<prod_id[i] << prod_id[j] <<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p, \"\"}," << endl;}
        
        // need to store information regarding the integrals involved in each product
        Bases_product product(prod_basis[prod_id[i]],prod_basis[prod_id[j]],prod_id[i],prod_id[j]);
        products_map[prod_id[i] + prod_id[j]] = product;
      
      }
    }
  
  }
  
  
  
  
  
  
  math_3 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_3.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_3.m ");
  if(options.verbose) system("./output/math_3.m");
  else system("./output/math_3.m  >/dev/null ");
#endif


}

void Calc_amplitudes::make_finite_amp(bool counter_terms)
{

  string prevb = "math_1.mx\"";
  string dimension = "D";
  ofstream math_4;
  math_4.open ("output/math_4.m");
  
  utils::print_math_header(math_4);
 // math_4<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n";
  
  print_finite_basis(reduced_basis,math_4);
  print_finite_basis(prod_basis,math_4);
  
  ////// new stuff below  //////
  
   // BASIS INTEGRAL COEFFICIENTS //
  format_coeff("D",reduced_basis,  reduced_basis_id, masses_input, id_input);
  //format_coeff_brackets(reduced_basis,  reduced_basis_id, masses_input, id_input);

  for (int i = 0; i < nbr;i++)
  {
    math_4 << "C" << reduced_basis_id[i] << " = " << reduced_basis[reduced_basis_id[i]].coefficient << ";" <<endl;
  }
  
  // PRODUCTS COEFFICIENTS //
  
  const char* file_integrals3_tmp = "output/output_products";
  string c_file_integrals3 = file_integrals3_tmp + blank + ext;
  const char *file_integrals3 = c_file_integrals3.c_str();
  vector<string> name_products,coeff_products_new;
  int temp_int;
  get_data(name_products, coeff_products_new, temp_int,file_integrals3, true);
  std::map <std::string, Bases > prod_map = products_container(prod_id);
  for (int i=0; i<temp_int; i++)
  {
    prod_map[name_products[i]].coefficient = coeff_products_new[i];
  }
  std::map <std::string, Bases > reduced_prod_map = remove_zeros(prod_map, extract_keys(prod_map));
  vector<string> reduced_prod_names = extract_keys(reduced_prod_map);

  format_coeff("D",reduced_prod_map,  reduced_prod_names, masses_input, id_input);
 // format_coeff_brackets(reduced_prod_map,  reduced_prod_names, masses_input, id_input);
  for (unsigned int i = 0; i < reduced_prod_names.size();i++)
  {
    math_4 << "C" << reduced_prod_names[i] << " = " << reduced_prod_map[reduced_prod_names[i]].coefficient << ";" <<endl;
  }
  
  //  SUMMATION
  
  math_4 << "  SEnFinite = 0 ";
  for (int i = 0; i<nbr;i++)
  {
     math_4  << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i];
  }
  for (unsigned int i = 0; i<reduced_prod_names.size();i++)
  {
    Bases temp_base;
    temp_base = reduced_prod_map[reduced_prod_names[i]];
    math_4 << " + "<< temp_base.e1 << " * " << temp_base.e2 << " * C" << reduced_prod_names[i];
  }
  math_4 << ";"<<endl;
  

  math_4<<"SEnFinite = SEnFinite /. D-> 4-2*epsilon;\n";
  math_4<<"Get[\"" << s_cwd <<"/output/remainder.mx\"]\n";
  math_4<<"SEn = SEnFinite + remainder;\n";
  
  if (counter_terms)
  {
    math_4<<"SEn = SEn * (1/epsilon)"<<endl;
  }
  
  math_4<<"SEn = Coefficient[SEn,epsilon,0]; \n";
  math_4<<"SEn = Simplify[SEn /. epsilon->0];\n"; // some integrals come through as D = 4-epsilon so fix these
  
  math_4<<"DumpSave[\""<<s_cwd<<"/output/math_2.mx\", SEn];\n";
  math_4<<"DumpSave[\""<<s_cwd<<"/output/math_ct.mx\", SEnFinite];\n";
  
  math_4.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_4.m ");
  if(options.verbose) system("./output/math_4.m");
  else system("./output/math_4.m  >/dev/null ");
#endif
}


void Calc_amplitudes::second_initial_trial(string prevb,string dimension)
{
  ofstream math_1;
  math_1.open ("output/math_1b.m");
  utils::print_math_header(math_1);
  math_1<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n";
  print_math_basis(full_basis,math_1,"SEn",dimension);
  math_1 << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;
  for (int i = 0; i < nb-1;i++)
  {
    math_1 << "{\""<<full_basis_id[i]<<" \", CForm[C"<<full_basis_id[i]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  math_1 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C"<<full_basis_id[nb-1]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  math_1 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_1.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_1b.m ");
  if (options.verbose) system("./output/math_1b.m");
  else system("./output/math_1b.m  >/dev/null");
#endif
  
}


void Calc_amplitudes::initial_trial(string dimension)
{
  ofstream math_1;
  math_1.open ("output/math_1.m");
  utils::print_math_header(math_1);
  utils::print_math_body(math_1,options,s_cwd,masses_input);
  math_1<<"Print[tfiamp0]\n"
  <<"SEn = Simplify[TarcerRecurse[tfiamp0] ];\n"
  <<"SEn = DiracSimplify[SEn];\n"
  <<"SEn = SEn /. DiracGamma[Momentum[p, D], D] -> p ;\n"
  <<"DumpSave[\""<<s_cwd<<"/output/math_1.mx\", SEn];\n"
  <<"Print[\"----------- The self energy is ---------- = \"]\n"
  <<"Print[SEn]\n"
  <<"Print[\"-------------------- = \"]"<< endl;
  
  print_math_basis(full_basis,math_1,"SEn",dimension);
  
  math_1 << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;
  for (int i = 0; i < nb-1;i++)
  {
    //math_1 << "{\""<<full_basis_id[i]<<" \", CForm[C"<<full_basis_id[i]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
    math_1 << "{\""<<full_basis_id[i]<<" \", C"<<full_basis_id[i]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p, \"\"}," << endl;
  }
  //math_1 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C"<<full_basis_id[nb-1]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  math_1 << "{\""<<full_basis_id[nb-1]<<" \", C"<<full_basis_id[nb-1]<<" /.Pair[Momentum[p],Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p, \"\"}" << endl;
  math_1 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_1.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_1.m ");
  if (options.verbose) system("./output/math_1.m");
  else system("./output/math_1.m  >/dev/null");
#endif
  
}




bool Calc_amplitudes::calc_diagram(Options options_in)
{
  
  bool success=0;
  bool sum_integrals=1;
  options = options_in;
  
  string diagram = options.diagram;
  string particle = options.particle;
  model = options.model;
  int loop_order = options.loop_order;
  
  
  cout << "calculating diagram " << diagram << " for particle " << particle << " in model ";
  cout << model << " at " << options.loop_order << "-loop order" << endl;
  if (options.counter_terms == true) { cout << "calculating counter term diagram"<< endl;}

  
  string particle_full = particle;
  
  particle =  part_name_simple(particle_full);
  // edit particle name to a safe string
  
  string tag="";
  if (options.counter_terms) {tag = particle + "_" + diagram + "_" + to_string(loop_order)+"c";}
  else {tag = particle + "_" + diagram + "_" + to_string(loop_order);}
  
   const char* file_masses_tmp = "models/";
  string c_file_masses = file_masses_tmp + model + "/masses" + ext;
  const char *file_masses = c_file_masses.c_str();
  int na;
  get_data(masses_input,id_input,na,file_masses);
  
  full_basis = set_bases(masses_input, id_input);
  full_basis_id = extract_keys(full_basis);
  nb = full_basis_id.size();
  
  
  // subroutine to generate and call Mathematica scripts
  initial_trial("D");
  compute_amp("math_1.mx\"","D");
  make_full_trial("math_1.mx\"","D",false);
  
  
  make_finite_amp(options.counter_terms);
  second_initial_trial("math_2.mx\"","4");
  compute_amp("math_2.mx\"","4");
  make_full_trial("math_2.mx\"","4",true);
  
  
  string remainder;
  success = check_done(remainder);
  
  
  // copy the Mathematica data file containing the full divergent amplitude
  const char *mx_ext = ".mx";
  const char* copy_tmp = "/output/math_data_";
  string c_copy = "cp output/math_ct.mx models/" + model +copy_tmp + tag + mx_ext;
  const char *copy_cmd = c_copy.c_str();
  system(copy_cmd);
  
  
  
  
  ReplaceAll(remainder,"Pair(Momentum(p),Momentum(p))", "Power(p,2)");
  
  // BASIS INTEGRAL COEFFICIENTS //
  
  const char* coeff_integrals_tmp = "/output/coeff_integrals_";
  string c_coeff_integrals = "models/" + model +coeff_integrals_tmp + tag + ext;
  const char *coeff_integrals = c_coeff_integrals.c_str();
  ofstream coeff_integrals_out;
  coeff_integrals_out.open (coeff_integrals);
  format_coeff("4",reduced_basis,  reduced_basis_id, masses_input, id_input);
  if (remainder != "0")
  {
    coeff_integrals_out << "  TSIL_COMPLEXCPP C0 = " << remainder << ";" <<endl;
  }
  for (int i = 0; i < nbr;i++)
  {
    coeff_integrals_out << "  TSIL_COMPLEXCPP C" << reduced_basis_id[i] << " = " << reduced_basis[reduced_basis_id[i]].coefficient << ";" <<endl;
  }
  coeff_integrals_out.close();
  
  
  // PRODUCTS COEFFICIENTS //
  
  const char* file_integrals3_tmp = "output/output_products";
  string c_file_integrals3 = file_integrals3_tmp + blank + ext;
  const char *file_integrals3 = c_file_integrals3.c_str();
  vector<string> name_products,coeff_products_new;
  int temp_int;
  get_data(name_products, coeff_products_new, temp_int,file_integrals3, true);
  std::map <std::string, Bases > prod_map = products_container(prod_id);
  for (int i=0; i<temp_int; i++)
  {
    prod_map[name_products[i]].coefficient = coeff_products_new[i];
  }
  std::map <std::string, Bases > reduced_prod_map = remove_zeros(prod_map, extract_keys(prod_map));
  vector<string> reduced_prod_names = extract_keys(reduced_prod_map);
  
  const char* coeff_products_tmp = "/output/coeff_products_";
  string c_coeff_products = "models/" + model +coeff_products_tmp + tag + ext;
  const char *coeff_products = c_coeff_products.c_str();
  
  ofstream coeff_products_out;
  coeff_products_out.open (coeff_products);
  format_coeff("4",reduced_prod_map,  reduced_prod_names, masses_input, id_input);
  for (unsigned int i = 0; i < reduced_prod_names.size();i++)
  {
    coeff_products_out << "  TSIL_COMPLEXCPP C" << reduced_prod_names[i] << " = " << reduced_prod_map[reduced_prod_names[i]].coefficient << ";" <<endl;
  }
  coeff_products_out.close();
  
  
  //  SUMMATION
  
  
  const char* summation_tmp = "/output/summation_"; // vector containing file names
  string c_summation = "models/" + model +summation_tmp + tag + ext;
  const char *summation = c_summation.c_str();
  
  ofstream summation_out;
  summation_out.open (summation);
  summation_out << "  return ";
  if (remainder != "0")
  {
    summation_out << " + C0 ";
  }
  
  
  for (int i = 0; i<nbr;i++)
  {
    if (sum_integrals != 0 ) summation_out  << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i];
  }
  for (unsigned int i = 0; i<reduced_prod_names.size();i++)
  {
    Bases temp_base;
    temp_base = reduced_prod_map[reduced_prod_names[i]];
    summation_out << " + "<< temp_base.e1 << " * " << temp_base.e2 << " * C" << reduced_prod_names[i];
  }
  
  if ((nbr==0) && (reduced_prod_names.size()==0) && (remainder == "0"))
  {
    summation_out << "0";
  }
  
  summation_out << ";"<<endl;
  summation_out.close();
  
  
  // BASIS INTEGRALS
  
  for (unsigned int i = 0; i < reduced_prod_names.size() ; i++)
  {
    reduced_basis_id.push_back(products_map[reduced_prod_names[i]].name_1);
    reduced_basis_id.push_back(products_map[reduced_prod_names[i]].name_2);
  }
  sort(reduced_basis_id.begin(),reduced_basis_id.end());
  reduced_basis_id.erase( unique( reduced_basis_id.begin(), reduced_basis_id.end() ), reduced_basis_id.end() );
  
  const char* basis_integrals_tmp = "/output/basis_integrals_"; // vector containing file names
  string c_basis_integrals = "models/" + model + basis_integrals_tmp + tag + ext;
  const char *basis_integrals = c_basis_integrals.c_str();
  
  ofstream basis_integral_out;
  basis_integral_out.open (basis_integrals);
  for (unsigned int i = 0; i<reduced_basis_id.size();i++)
  {
    basis_integral_out << reduced_basis_id[i] << endl;
  }
  basis_integral_out.close();
  
  if (success) {update_avail_diagrams(options);}
  return success;
}


void draw_all_diagrams(Options options)
{
  string s_cwd(getcwd(NULL,0));
  string particle = options.particle;
  string model = options.model;
  
  ofstream myfile;
  myfile.open ("output/make_figures.m");
  
  string type="";
  utils::print_math_header(myfile);
  if (options.counter_terms){myfile<<"t12 = CreateCTTopologies["<<options.loop_order<<", 1 -> 1, ExcludeTopologies -> Internal];"<<endl;
    type = to_string(options.loop_order) +"c";}
  else {myfile<<"t12 = CreateTopologies["<<options.loop_order<<", 1 -> 1, ExcludeTopologies -> Internal];"<<endl;
    type = to_string(options.loop_order);}
  myfile<<"alldiags = InsertFields[t12, {"<<particle<<"} -> {"<<particle<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<s_cwd<<"/models/"<<model<<"/"<<model<<"\"];\n"
  <<"Export[\""<<s_cwd<<"/models/"<<options.model<<"/FA_diagrams/diagrams_"<<particle<< "_" << type <<".pdf\",Paint[alldiags,Numbering->Simple]];\n"  // print the FA diagram to pdf in local directory
  <<endl;
  
#ifdef RUN_ALL
  system("chmod +x output/make_figures.m ");
  if (options.verbose) system("./output/make_figures.m");
  else system("./output/make_figures.m  >/dev/null");
#endif
}




void draw_diagrams(vector<std::string> particles, vector<std::string> diagrams, int nd,string model)
{
  
  string s_cwd(getcwd(NULL,0));
  ofstream myfile;
  myfile.open ("output/make_figures.m");
  
  Options options;
  
  utils::print_math_header(myfile);
  
  myfile <<"t12 = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> Internal];\n"
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
      }
    }
    myfile <<"]\n";
    myfile <<"Export[\""<<s_cwd<<"/FA_diagrams/subset_diagrams_"<<particle_name_tmp<<".pdf\",Paint[subdiags"<<particle_name_tmp<<"]];\n"  // print the FA diagram to pdf in local directory
    <<endl;
  }
  
#ifdef RUN_ALL
  system("chmod +x output/make_figures.m ");
  if (options.verbose) system("./output/make_figures.m");
  else system("./output/make_figures.m  >/dev/null");
#endif
  
  
  
}


void Calc_amplitudes::generate_figures(Options options_in)
{
  
  vector<std::string> particles, diagrams;
  
  options = options_in;
  
  
  if (options.input_list!="")
  {
    std::string particles_in [1000];
    std::string diagrams_in [1000]; int i=0;
    std::ifstream input(options.input_list);
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
  else
  {
    draw_all_diagrams(options);
  }
  
  
}
