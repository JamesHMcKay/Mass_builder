/*
 Mass Builder
 
 James McKay
 Aug 2016 - June 2017
 
 --- compute_amp.cpp ---
 
 This file contains the main routine for calling FeynArts, FeynCalc and TARCER,
 followed by sorting this output and storing it for later use in generated
 TSIL interface.
 
 This file also includes functions for sending FeynArts diagrams directory to pdf in the folder FAdiagrams
 and the functions for computing the tree-level counter-term coupling
 
 */

#include "compute_amp.hpp"

//#define DEBUG
using namespace std;
using namespace utils;
using namespace templates;


// main routine for computing an amplitude for a given diagram
// returns true for success and false if an error is encountered
bool Compute_amp::calc_diagram()
{
  bool success = 0;
  
  // print diagram details (number, particle, ...) to terminal
  print_diagram_info(options);
  
  string tag = make_tag(options);
  
  // get list of masses and identifier strings for this model
  const char* file_masses_tmp = "models/";
  string c_file_masses = file_masses_tmp + options.model + "/masses" + ext;
  const char *file_masses = c_file_masses.c_str();
  int na;
  get_data(masses_input,id_input,na,file_masses);
  
  // create WSTP link and compute amplitude
  create_wstp_link();
  load_libraries();
  
  // create string containing input to Mathematica
  std::string input;
  
  
  
  // check if this diagram has already been computed (need to check available diagrams list and the existence of the math file)
  
  if (check_if_available(options) && !options.force)
  {
    cout << "using previously computed amplitude " << endl;
    templates::print_math_header(input);
    utils::get_saved_amplitude(input,options,masses_input[0]);
    send_to_math(input);
  }
  
  else
  {
    cout << "computing amplitude " << endl;
    templates::print_math_header(input);
    utils::print_math_body(input,options,get_cwd(),masses_input);
    // send input to Mathematica
    send_to_math(input);
    
    // send command to apply TARCER recurse function
    utils::print_tarcer_recurse(input);
    send_to_math(input);
    
    // expand basis integrals into finite plus divergent pieces (using MassBuilder. package)
    input += "SelfEnergyFinite = expandBasisIntegrals[SE, masses, massesExpand, MassBuilderA,";
    input += "MassBuilderB, MassBuilderJ, MassBuilderK, MassBuilderT, MassBuilderV, MassBuilderF];";
    
    if (options.counter_terms && options.loop_order == 1)
    {
      // add 1/epsilon^2 order counter-terms to the tree-level counter-term amplitude since FeynArts doesn't include these
      if (options.fermion)
      {
        input += "SelfEnergyFinite = addHigherOrderDivergencesFermion[SelfEnergyFinite];";
      }
      else
      {
        input += "SelfEnergyFinite = addHigherOrderDivergences[SelfEnergyFinite];";
      }
    }
    
    if (options.counter_terms)
    {
      input += "SelfEnergyFinite = SelfEnergyFinite/MassBuilderEpsilon;";
    }
    
    // save full expanded amplitude to a Mathematica data file
    input += "DumpSave[\"" + get_cwd() + "/models/" + options.model + "/output/math_data_" + tag + ".mx\", SelfEnergyFinite];";
    send_to_math(input);
  }
  
  get_finite_amplitude(input,options);
  
  send_to_math(input);
  
  // extract non-zero coefficients
  // send List of all possible basis integrals and get list back of coefficients
  
  // create std::map<std::string, Bases> containing all Bases objects
  full_basis = set_bases(masses_input, id_input);
  // obtain a corresponding vector<string> of unique Bases identifiers
  full_basis_id = extract_keys(full_basis);
  nb = full_basis_id.size();
  
  // generate a .m file to call within Mathematica (since it is a very long list)
  ofstream math_1;
  math_1.open (add_mpi_ext("output/math_1", options.mpi_process, "m"));
  print_math_basis(full_basis,math_1,"SelfEnergyFinite","4");
  math_1 << "Export[\""<<get_cwd()<<"/output/output_"<< options.mpi_process << ".txt\", {" << endl;
  
  for (int i = 0; i < nb-1;i++)
  {
    math_1 << "{\""<<full_basis_id[i]<<" \", CForm[C"<<full_basis_id[i]<<" /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  math_1 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C"<<full_basis_id[nb-1]<<"  /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  math_1 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_1.close();
  
  
  // call the .m file in Mathematica
  input+= "AppendTo[$Path, \"" + get_cwd() + "/output/\"];";
  string filename = add_mpi_ext("output/math_1", options.mpi_process, "m");
  input += "<< " + filename + ";";
  
  send_to_math(input);
  
  // now attempt to construct amplitude //
  
  // now dealing a with finite (4 dimensional) amplitude
  string dimension = "4";
  vector<std::string> output_string, coeff_new;
  int temp_int = 0;
  
  // get the output from Mathematica (list of coefficients and bases object identifiers)
  const char* file_integrals2_tmp = "output/output_";
  string c_file_integrals2 = file_integrals2_tmp + std::to_string(options.mpi_process) + ext;
  const char *file_integrals2 = c_file_integrals2.c_str();
  get_data(output_string, coeff_new, temp_int,file_integrals2, true);
  
  // assign the coefficients computed by Mathematica to each Bases object
  for (int i=0; i<temp_int; i++)
  {
    full_basis[output_string[i]].coefficient = coeff_new[i];
  }
  
  // discard all the Bases objects that have a zero coefficient
  reduced_basis  = remove_zeros(full_basis, full_basis_id);
  reduced_basis_id = extract_keys(reduced_basis);
  nbr = reduced_basis_id.size();
  
  //  done updating basis map //
  
  // we construct the trial amplitude using the reduced list of required integrals
  
  // the trail amplitude construction is done in a .m file which we will call in Mathematica
  ofstream math_2;
  math_2.open (add_mpi_ext("output/math_2", options.mpi_process, "m"));
  
  // print the bases to the .m file
  print_math_basis(reduced_basis,math_2, "SelfEnergyFinite","4");
  
  // construct the amplitude as coefficients multiplied by basis integrals
  math_2 << "SelfEnergyTrial = 0 ";
  for (int i = 0; i<nbr;i++)
  {
    math_2 << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i]; // use the terms we know are non-zero
  }
  math_2 << ";"<<endl;
  
  // evaluate the difference between the trial and actual amplitude
  math_2 << "difference = Simplify[SelfEnergyFinite-SelfEnergyTrial];"<<endl;

  // save the difference to a text file for reference
  math_2 << "Export[\""<<get_cwd()<<"/output/result_"<< options.mpi_process << ".txt\", difference];" << endl;
  
  // evaluate the coefficients for the difference (see what basis integrals appear in the difference contains)
  // and write this to file
  print_math_basis(full_basis,math_2, "difference",dimension);
  math_2 << "Export[\""<<get_cwd()<<"/output/output2_"<< options.mpi_process << ".txt\", {" << endl;
  for (int i = 0; i < nb-1;i++)
  {
    math_2 << "{\""<<full_basis_id[i]<<" \", CForm[C"<<full_basis_id[i]<<" + C" <<full_basis_id[i] <<"2  /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  math_2 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C"<<full_basis_id[nb-1]<<" + C"<< full_basis_id[nb-1] <<"2  /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
  math_2 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_2.close();
  
  
  // ask Mathematica to evaluate math_2.m
  filename = add_mpi_ext("output/math_2", options.mpi_process, "m");
  input += "<< " + filename + ";";
  send_to_math(input);
  
  // read in the result of the above check
  temp_int = 0;
  const char* file_integrals4_tmp = "output/output2_";
  string c_file_integrals4 = file_integrals4_tmp + std::to_string(options.mpi_process) + ext;
  const char *file_integrals4 = c_file_integrals4.c_str();
  vector<std::string> coeff_new_prod;
  prod_basis.clear();
  prod_id.clear();
  get_data(prod_id, coeff_new_prod, temp_int,file_integrals4,true);
  
  // create a fresh set of Bases objects (from the original full set)
  // but without any F objects, since we know these won't appear in the
  // difference (the difference only contains products of basis integrals)
  prod_basis = remove_type_F(full_basis,full_basis_id);
  
  // check that we got some input back
  if (temp_int == 0)
  {
    cout << "ERROR LIST OF BASIS INTEGRALS IS EMPTY" << endl;
  }
  else
  {
    // assign the coefficients to the Bases objects
    for (int i=0; i<temp_int; i++)
    {
      prod_basis[prod_id[i]].coefficient = coeff_new_prod[i];
    }
    // remove the Bases objects that have a zero coefficient
    prod_basis = remove_zeros(prod_basis, prod_id);
    prod_id = extract_keys(prod_basis);
    
    // this should be redundant?
    prod_basis = remove_type_F(prod_basis, prod_id);
    prod_id = extract_keys(prod_basis);
  }
  
  // if no products were found but amplitude not fully constructed
  // then set products basis to be all, this should not occur
  // return error here? (todo: check if this ever happens)
  if (prod_basis.size()==0 && !check_done_quiet(options.mpi_process))
  {
    prod_basis = reduced_basis;
    prod_id = extract_keys(prod_basis);
    prod_basis = remove_type_F(prod_basis, prod_id);
    prod_id = extract_keys(prod_basis);
  }
  
  np = prod_id.size();
  
  
  
  // generate final .m file to construct the amplitude
  ofstream math_3;
  math_3.open (add_mpi_ext("output/math_3", options.mpi_process, "m"));
  
  print_math_basis(reduced_basis,math_3,"SelfEnergyFinite","4");
  print_math_basis(prod_basis,math_3,"SelfEnergyFinite","4");
  print_math_products(prod_basis,math_3,"SelfEnergyFinite","4");
  
  math_3 << "SelfEnergyTrial = 0 ";
  // consider basis integrals that are not products
  for (int i = 0; i<nbr;i++)
  {
    math_3 << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i];
  }
  // consider the products of basis integrals
  for (int i = 0; i<np ; i++)
  {
    for (int j = 0; j<np ; j++)
    {
      math_3 << " + "<< prod_id[i] << prod_id[j] << " * C"<<prod_id[i]<< prod_id[j];
    }
  }
  math_3 << ";"<<endl;
  
  // evaluate the difference between the constructed amplitude and the original amplitude
  math_3 << "diff = Simplify[SelfEnergyFinite-SelfEnergyTrial];\n";

  // save the result to a text file
  math_3 << "Export[\""<<get_cwd()<<"/output/result_"<< options.mpi_process << ".txt\", CForm[diff/. DiracGamma[Momentum[p]] -> p] ]" << endl;
  math_3<< "remainder = diff;"<<endl;
  // save the remainder to a Mathematica data file
  math_3<< "DumpSave[\""<<get_cwd()<<"/output/remainder_"<< options.mpi_process << ".mx\", remainder];"<<endl;
  math_3 << "Export[\""<<get_cwd()<<"/output/output_products_"<< options.mpi_process << ".txt\", {" << endl;
  
  products_map.clear();
  
  bool cform = true;
  // write out the final coefficients in C++ form
  for (int i = 0; i<np;i++)
  {
    for (int j = 0; j<np;j++)
    {
      if (!cform)
      {
        if ((i==np-1) && (j==np-1)){ math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",C"<<prod_id[i] << prod_id[j] <<"  /. DiracGamma[Momentum[p]] -> p, \"\"}" << endl;}
        else {math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",C"<<prod_id[i] << prod_id[j] <<" /. DiracGamma[Momentum[p]] -> p, \"\"}," << endl;}
      }
      else
      {
        if ((i==np-1) && (j==np-1)){ math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C"<<prod_id[i] << prod_id[j] <<"  /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;}
        else {math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C"<<prod_id[i] << prod_id[j] <<" /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;}
      }
      Bases_product product(prod_basis[prod_id[i]],prod_basis[prod_id[j]],prod_id[i],prod_id[j]);
      products_map[prod_id[i] + prod_id[j]] = product;
    }
  }

  math_3 << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  math_3.close();
  
  filename = add_mpi_ext("output/math_3", options.mpi_process, "m");
  input += "<< " + filename + ";";
  
  send_to_math(input);
  
  input += "Quit[]";
  send_to_math(input);
  
  WSClose(link);
  cout << "WSTP link closed successfully" << endl;
  
  
  ///// organise output data////////
  
  // check if the remainder contains any basis integrals
  string remainder;
  success = check_done(remainder,options.mpi_process);
  

  
  ReplaceAll(remainder,"Pair(Momentum(p),Momentum(p))", "(-Power(p,2))");
  
  // BASIS INTEGRAL COEFFICIENTS //
  
  const char* coeff_integrals_tmp = "/output/coeff_integrals_";
  string c_coeff_integrals = "models/" + options.model +coeff_integrals_tmp + tag + ext;
  const char *coeff_integrals = c_coeff_integrals.c_str();
  ofstream coeff_integrals_out;
  coeff_integrals_out.open (coeff_integrals);
  format_coeff("4",reduced_basis,  reduced_basis_id, masses_input, id_input);
  if (remainder != "0")
  {
    format_coeff(remainder,masses_input,id_input);
    coeff_integrals_out << "  TSIL_COMPLEXCPP C0 = " << remainder << ";" <<endl;
  }
  for (int i = 0; i < nbr;i++)
  {
    coeff_integrals_out << "  TSIL_COMPLEXCPP C" << reduced_basis_id[i] << " = " << reduced_basis[reduced_basis_id[i]].coefficient << ";" <<endl;
  }
  coeff_integrals_out.close();
  
  
  // PRODUCTS COEFFICIENTS //
  
  const char* file_integrals3_tmp = "output/output_products_";
  string c_file_integrals3 = file_integrals3_tmp + std::to_string(options.mpi_process) + ext;
  const char *file_integrals3 = c_file_integrals3.c_str();
  vector<string> name_products,coeff_products_new;
  
  get_data(name_products, coeff_products_new, temp_int,file_integrals3, true);
  std::map <std::string, Bases > prod_map = products_container(prod_id);
  for (int i=0; i<temp_int; i++)
  {
    prod_map[name_products[i]].coefficient = coeff_products_new[i];
  }
  std::map <std::string, Bases > reduced_prod_map = remove_zeros(prod_map, extract_keys(prod_map));
  vector<string> reduced_prod_names = extract_keys(reduced_prod_map);
  
  const char* coeff_products_tmp = "/output/coeff_products_";
  string c_coeff_products = "models/" + options.model +coeff_products_tmp + tag + ext;
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
  string c_summation = "models/" + options.model +summation_tmp + tag + ext;
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
    summation_out  << " + "<<reduced_basis_id[i]<< " * C"<<reduced_basis_id[i];
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
  string c_basis_integrals = "models/" + options.model + basis_integrals_tmp + tag + ext;
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




std::string draw_all_diagrams(std::string &input, Options options)
{
  string particle_1 = options.particle_1;
  string particle_2 = options.particle_2;
  string tag = particle_1;
  
  if (particle_1 != particle_2)
  {
    tag = particle_1 + "_" + particle_2;
  }
  
  int n_final_states = options.n_final_states;
  
  string type="";
  templates::print_math_header(input);
  
  utils::assign_FCGV(input,options);
  utils::assign_variables(input,options);
  
  if (options.counter_terms)
  {
    input+= "t12 = CreateCTTopologies[" + to_string(options.loop_order) + ", 1 -> 1, ExcludeTopologies -> Internal];";
    type = to_string(options.loop_order) +"c";
  }
  else
  {
    input+= "t12 = CreateTopologies[" + to_string(options.loop_order) + ", 1 -> " +  to_string(n_final_states)  +  ", ExcludeTopologies -> Internal];";
    type = to_string(options.loop_order);
  }
  
  if (options.use_lorentz)
  {
    input+= "alldiags = InsertFields[t12, {" + particle_1 + "} -> {" + particle_2 + "},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions -> {"  +  options.restrictions  +  "},Model -> \"" + get_cwd() + "/models/" + options.model + "/" + options.model + "\"];\n";
  }
  else
  {
    input+= "alldiags = InsertFields[t12, {" + particle_1 + "} -> {" + particle_2 + "},InsertionLevel -> {Particles}, GenericModel -> \"" + get_cwd() + "/models/" + options.model + "/" + options.model + "\",Restrictions -> {"  +  options.restrictions  +  "},Model -> \"" + get_cwd() + "/models/" + options.model + "/" + options.model + "\"];";
  }
  
  input+= "Export[\"" + get_cwd() + "/models/" + options.model + "/FA_diagrams/diagrams_"+ tag + "_" + type + ".pdf\",Paint[alldiags,Numbering->Simple]];";  // print the FA diagram to pdf in local directory
  
  return input;
}


void Compute_amp::generate_figures()
{
  vector<std::string> particles, diagrams;
  
  std::string input;
  draw_all_diagrams(input,options);
  create_wstp_link();
  load_libraries();
  send_to_math(input);
  input += "Quit[]";
  send_to_math(input);
  
  WSClose(link);
  cout << "WSTP link closed successfully" << endl;
  
  }





void Compute_amp::solve_1loop(std::string particle,vector<std::string> diagram)
{
  
  
  int nd = diagram.size();
  
  string dimension = "D";
  
  create_wstp_link();
  load_libraries();
  WSNewPacket(link);
  
  std::string input;
  

  
  input += "  SEtotal = 0 ;";
  
  for (int i=0;i<nd;i++)
  {
    
    input += "Get[\"" + get_cwd() + "/models/" + options.model + "/output/math_data_" + particle + "_" + diagram[i] + "_1" + ".mx\"];";
    
    input += "SEtotal = SEtotal + SelfEnergyFinite;";
  }
  
  input += "SEtotal = SEtotal /. D -> 4 - 2*MassBuilderEpsilon;";
  input += "SE = Coefficient[SEtotal,MassBuilderEpsilon,-1];";
  input += "SE = Simplify[SE /. MassBuilderEpsilon->0];";
  
  
  // check for higher orders in 1/epsilon
  input += "SEhot = Coefficient[SEtotal,MassBuilderEpsilon,-2] + Coefficient[SEtotal,MassBuilderEpsilon,-3];";
  
  input += "ToString[SEhot,CForm];";
  
  send_to_math(input);
  
  const char* high_order_terms;
  
  if(!WSGetString((WSLINK)pHandle, &high_order_terms))
  {
    cout << "Error getting string from WSTP" << endl;
  }
  
  
  
  cout << "Higher order divergences in (1/epsilon) are = " << high_order_terms << endl;
  
  input += "ct = FullSimplify[-SE*Pi^2/.Pair[Momentum[p], Momentum[p]]->p^2];";
  
  input += "ToString[ct,CForm]";
  
  send_to_math(input);
  
  const char* counter_term;
  
  if(!WSGetString((WSLINK)pHandle, &counter_term))
  {
    cout << "Error getting string from WSTP" << endl;
  }
  
  cout << "Counter-term coupling = " << counter_term << endl;
  
  input += "Quit[]";
  send_to_math(input);
  
  WSClose(link);
  cout << "WSTP link closed successfully" << endl;
  

}





void Compute_amp::calc_counter_terms()
{
  // need to read in the list of available diagrams and then select the 1-loop ones for
  // adding up here
  
  // read in available diagrams
  
  const char *ext = ".txt";
  const char* file_diagrams_tmp = "models/";
  string c_file_diagrams = file_diagrams_tmp + options.model + "/output/avail_diagrams_" + ext;
  const char *file_diagrams = c_file_diagrams.c_str();
  
  vector<std::string> tags;
  vector<std::string> particle_names,levels;
  string level;
  int nd; // number of diagrams
  
  cout << "input list = " << file_diagrams << endl;
  
  get_data(particle_names,tags,levels, nd,file_diagrams);
  
  // one-loop corrections
  vector<std::string> tags_1;
  vector<std::string> particle_names_1,levels_1;
  
  for (int i=0;i<nd;i++)
  {
    if (levels[i]=="1" && particle_names[i]==options.particle)
    {
      particle_names_1.push_back(particle_names[i]);
      tags_1.push_back(tags[i]);
      levels_1.push_back(levels[i]);
    }
    
  }
  
  string particle_simple = part_name_simple(options.particle_1,options.particle_2);
  
  solve_1loop(particle_simple,tags_1);
  
}

