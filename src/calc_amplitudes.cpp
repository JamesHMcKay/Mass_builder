/*
 Mass Builder 
 
 James McKay
 Aug - Sep 2016
 
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


bool Calc_amplitudes::calc_diagram(Options options)
{
  bool success=0;
  bool sum_integrals=1;
  
  string diagram = options.diagram;
  string particle = options.particle;
  string model = options.model;
  int loop_order = options.loop_order;
  
  cout << "calculating diagram " << diagram << " for particle " << particle << " in model ";
  cout << model << " at " << options.loop_order << "-loop order" << endl;
  if (options.counter_terms == true) { cout << "calculating counter term diagram"<< endl;}
  
  const char *ext = ".txt";
  string underscore = "_";
  string blank = "";
  
  string particle_full = particle;
  
  particle =  part_name_simple(particle_full);
  // edit particle name to a safe string
  
  string tag="";
  if (options.counter_terms) {tag = particle + "_" + diagram + "_" + to_string(loop_order)+"c";}
  else {tag = particle + "_" + diagram + "_" + to_string(loop_order);}
  
  int n = 0;
  vector<std::string> output0, output1,output2;
  vector<double> dup;
  vector<std::string> A;
  
  const char* file_masses_tmp = "models/";
  string c_file_masses = file_masses_tmp + model + "/masses" + ext;
  const char *file_masses = c_file_masses.c_str();
  
  
  get_data(A, n, file_masses);
  
  string s_cwd(getcwd(NULL,0));
  
  output0.resize(pow(n,5));
  
  vector<string> masses_input,id_input;
  int na;
  get_data(masses_input,id_input,na,file_masses);
  std::map<std::string, Bases> base_map = set_bases(masses_input, id_input);
  
  vector<string> bases_names = extract_keys(base_map);
  
  ofstream myfile;
  myfile.open ("output/math_1.m");
  
  utils::print_math_header(myfile);
  utils::print_math_body(myfile,options,s_cwd,masses_input);
  
  myfile<<"Print[tfiamp0]\n"
  <<"SEn = Simplify[TarcerRecurse[tfiamp0] ];\n"
  <<"DumpSave[\""<<s_cwd<<"/output/math_1.mx\", SEn];\n"
  <<"Print[\"----------- The self energy is ---------- = \"]\n"
  <<"Print[SEn]\n"
  <<"Print[\"-------------------- = \"]"<< endl;
  
  print_math_basis(base_map,myfile, "SEn");
  
  myfile << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;
  for (unsigned int i = 0; i < bases_names.size()-1;i++)
  {
    myfile << "{\""<<bases_names[i]<<" \", CForm[C"<<bases_names[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  myfile << "{\""<<bases_names[bases_names.size()-1]<<" \", CForm[C"<<bases_names[bases_names.size()-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  myfile << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  myfile.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_1.m ");
  if (options.verbose) system("./output/math_1.m");
  else system("./output/math_1.m  >/dev/null");
#endif
  
  
  ///  part 1 complete ////
  
  ///  part 2 below ////
  
  
  
  
  vector<std::string> output_string, coeff_new;
  int temp_int = 0;
  
  const char* file_integrals2_tmp = "output/output";
  string c_file_integrals2 = file_integrals2_tmp + blank + ext;
  const char *file_integrals2 = c_file_integrals2.c_str();
  
  get_data(output_string, coeff_new, temp_int,file_integrals2, true);
  for (int i=0; i<temp_int; i++)
  {
    base_map[output_string[i]].coefficient = coeff_new[i]; // set coefficients from first math run
  }
  
  std::map <std::string, Bases > reduced_base_map = remove_zeros(base_map, bases_names); // remove integrals with zero coefficients
  vector<string> reduced_bases_names = extract_keys(reduced_base_map);
  
  string prevb = "math_1.mx\"";
  ofstream math_2;
  math_2.open ("output/math_2.m");
  
  utils::print_math_header(math_2);
  math_2<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n"
  <<"Print[SEn]"<< endl;
  
  print_math_basis(reduced_base_map,math_2, "SEn");
  
  math_2 << "SEnTrial = 0 ";
  
  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
    
    math_2 << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
    
  }
  math_2 << ";"<<endl;
  
  math_2 << "diff = Simplify[SEn-SEnTrial]"<<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  math_2 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  math_2 << "Print[\" --------------------------------------- \"]" <<endl;
  
  math_2 << "Print[diff]"<<endl;
  
  math_2 << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;
  
  print_math_basis(reduced_base_map,math_2, "diff");
  
  math_2 << "Export[\""<<s_cwd<<"/output/output2.txt\", {" << endl;
  
  
  if (reduced_bases_names.size()!=0)
  {
    for (unsigned int i=0;i<reduced_bases_names.size()-1;i++)
    {
      math_2 << "{\""<<reduced_bases_names[i]<<" \", CForm[C"<<reduced_bases_names[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
    }
    
    math_2 << "{\""<<reduced_bases_names[reduced_bases_names.size()-1]<<" \", CForm[C"<<reduced_bases_names[reduced_bases_names.size()-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
    math_2 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  }
  
  if (reduced_bases_names.size()==0)
  {
    math_2 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  }
  math_2 << "Print[\"Completed\"]"<<endl;
  
  math_2.close();
  
  
  
#ifdef RUN_ALL
  system("chmod +x output/math_2.m ");
  if(options.verbose) system("./output/math_2.m ");
  else system("./output/math_2.m  >/dev/null ");
#endif
  
  
  // part 2 complete -- first attempt at constructing amplitude //
  
  // part 3 below, deal with products ///
  
  
  
  n = 0;
  vector<std::string> output_string_prod, coeff_new_prod;
  std::map <std::string, Bases > base_map_prod;
  std::map <std::string, Bases > reduced_base_map_copy=reduced_base_map; // copy the reduced map of basis integrals (should this be the original base map?)
  vector<string> bases_names_prod;
  temp_int = 0;
  
  const char* file_integrals4_tmp = "output/output2";
  string c_file_integrals4 = file_integrals4_tmp + blank + ext;
  const char *file_integrals4 = c_file_integrals4.c_str();
  
  get_data(output_string_prod, coeff_new_prod, temp_int,file_integrals4,true);
  
  if (temp_int == 0)
  {
    for (unsigned int i=0; i<bases_names.size(); i++)
    {
      if (base_map[bases_names[i]].type != "F")
      {
        base_map_prod[bases_names[i]] = base_map[bases_names[i]];
        bases_names_prod.push_back(bases_names[i]);
      }
    }
  }
  else
  {
    for (int i=0; i<temp_int; i++)
    {
      reduced_base_map_copy[output_string_prod[i]].coefficient = coeff_new_prod[i];
    }
    base_map_prod = remove_zeros(reduced_base_map_copy, reduced_bases_names);
    bases_names_prod = extract_keys(base_map_prod);
  }
  
  ofstream math_3;
  math_3.open ("output/math_3.m");
  
  prevb = "math_1.mx\"";
  utils::print_math_header(math_3);
  math_3<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n"
  <<"Print[SEn]"<< endl;
  
  
  if (temp_int !=0) {print_math_basis(reduced_base_map,math_3,"SEn");}
  else {print_math_basis(base_map,math_3,"SEn");}
  
  print_math_products(base_map_prod,math_3,"SEn");
  
  math_3 << "SEnTrial = ";
  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
    math_3 << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
  }
  
  for (unsigned int i = 0; i<bases_names_prod.size();i++)
  {
    for (unsigned int j = 0; j<bases_names_prod.size();j++)
    {
      math_3 << " + "<< bases_names_prod[i] << bases_names_prod[j] << " * C"<<bases_names_prod[i]<< bases_names_prod[j];
    }
  }
  math_3 << ";"<<endl;
  
  math_3 << "diff = Simplify[SEn-SEnTrial]"<<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  math_3 << "Print[\" --------------------------------------- \"]" <<endl;
  
  math_3 << "Print[diff]"<<endl;
  math_3 << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;
  
  // print out coefficients of products
  int nn = bases_names_prod.size();
  math_3 << "Export[\""<<s_cwd<<"/output/output_products.txt\", {" << endl;
  
  std::map <std::string, Bases_product > products_map;
  
  for (int i = 0; i<nn;i++)
  {
    for (int j = 0; j<nn;j++)
    {
      
      if ((i==nn-1) && (j==nn-1)){ math_3 << "{\""<<bases_names_prod[i] << bases_names_prod[j] <<" \",CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;}
      else {math_3 << "{\""<<bases_names_prod[i] << bases_names_prod[j] <<" \",CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
      
      // need to store information regarding the integrals involved in each product
      
      Bases_product product(base_map_prod[bases_names_prod[i]],base_map_prod[bases_names_prod[j]],bases_names_prod[i],bases_names_prod[j]);
      
      products_map[bases_names_prod[i] + bases_names_prod[j]] = product;
      
    }
  }
  
  math_3 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  
  
  math_3.close();
  
#ifdef RUN_ALL
  system("chmod +x output/math_3.m ");
  if(options.verbose) system("./output/math_3.m");
  else system("./output/math_3.m  >/dev/null ");
#endif
  
  success = check_done();
  
  
  ofstream diagram_list;
  diagram_list.open( "output/diagrams.txt", ios::out | ios::app );
  diagram_list << particle << " " << diagram << endl;
  diagram_list.close();
  
  const char* coeff_integrals_tmp = "/output/coeff_integrals_"; // vector containing file names
  string c_coeff_integrals = "models/" + model +coeff_integrals_tmp + tag + ext;
  const char *coeff_integrals = c_coeff_integrals.c_str();
  
  ofstream coeff_integrals_out;
  
  coeff_integrals_out.open (coeff_integrals);
  
  format_coeff(reduced_base_map,  reduced_bases_names, masses_input, id_input);
  
  for (unsigned int i = 0; i < reduced_bases_names.size();i++)
  {
    coeff_integrals_out << "TSIL_COMPLEXCPP C" << reduced_bases_names[i] << " = " << reduced_base_map[reduced_bases_names[i]].coefficient << ";" <<endl;
  }
  
  // list of product coefficients in form "TSIL_COMPLEX name = coefficient"
  // create simple containers for products
  
  std::map <std::string, Bases > prod_map = products_container(bases_names_prod);
  
  
  const char* file_integrals3_tmp = "output/output_products";
  string c_file_integrals3 = file_integrals3_tmp + blank + ext;
  const char *file_integrals3 = c_file_integrals3.c_str();
  vector<string> name_products,coeff_products_new;
  
  get_data(name_products, coeff_products_new, temp_int,file_integrals3, true);
  
  for (int i=0; i<temp_int; i++)
  {
    prod_map[name_products[i]].coefficient = coeff_products_new[i];
  }
  
  std::map <std::string, Bases > reduced_prod_map = remove_zeros(prod_map, extract_keys(prod_map));
  vector<string> reduced_prod_names = extract_keys(reduced_prod_map);
  
 
  
  
  
  const char* coeff_products_tmp = "/output/coeff_products_"; // vector containing file names
  string c_coeff_products = "models/" + model +coeff_products_tmp + tag + ext;
  const char *coeff_products = c_coeff_products.c_str();
  
  ofstream coeff_products_out;
  
  coeff_products_out.open (coeff_products);
  
  format_coeff(reduced_prod_map,  reduced_prod_names, masses_input, id_input);
  
  for (unsigned int i = 0; i < reduced_prod_names.size();i++)
  {
    coeff_products_out << "TSIL_COMPLEXCPP C" << reduced_prod_names[i] << " = " << reduced_prod_map[reduced_prod_names[i]].coefficient << ";" <<endl;
  }
  
  
  //  SUMMATION
  
  
  const char* summation_tmp = "/output/summation_"; // vector containing file names
  string c_summation = "models/" + model +summation_tmp + tag + ext;
  const char *summation = c_summation.c_str();
  
  ofstream summation_out;
  
  summation_out.open (summation);
  
  summation_out << "return ";
  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
    if (sum_integrals != 0 ) summation_out  << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
  }
  for (unsigned int i = 0; i<reduced_prod_names.size();i++)
  {
    
    Bases temp_base;
    temp_base = reduced_prod_map[reduced_prod_names[i]];
    summation_out << " + "<< temp_base.e1 << " * " << temp_base.e2 << " * C" << reduced_prod_names[i];
    
  }
  
  summation_out << ";"<<endl;
  summation_out.close();
  
  
  // Basis integrals list
  // get the integrals associated with the non-zero product coefficients
  
  for (unsigned int i = 0; i < reduced_prod_names.size() ; i++)
  {
    reduced_bases_names.push_back(products_map[reduced_prod_names[i]].name_1);
    reduced_bases_names.push_back(products_map[reduced_prod_names[i]].name_2);
  }
  
  sort(reduced_bases_names.begin(),reduced_bases_names.end());
  reduced_bases_names.erase( unique( reduced_bases_names.begin(), reduced_bases_names.end() ), reduced_bases_names.end() );
  
  
  const char* basis_integrals_tmp = "/output/basis_integrals_"; // vector containing file names
  string c_basis_integrals = "models/" + model + basis_integrals_tmp + tag + ext;
  const char *basis_integrals = c_basis_integrals.c_str();
  
  
  ofstream output_1;
  output_1.open (basis_integrals);
  
  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
    output_1 << reduced_bases_names[i] << endl;
  }
  
  output_1.close();

  
  
  
  update_avail_diagrams(options);
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


void Calc_amplitudes::generate_figures(Options options)
{
  
  vector<std::string> particles, diagrams;
  string model;
  
  
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
