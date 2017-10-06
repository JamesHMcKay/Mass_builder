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
  
	open_log_files();
	
  // print diagram details (number, particle, ...) to terminal
  print_diagram_info(options);
  
  string tag = make_tag(options);
  
  // get list of masses and identifier strings for this model
  const char* file_masses_tmp = "models/";
  string c_file_masses = file_masses_tmp + options.model + "/masses" + ext;
  const char *file_masses = c_file_masses.c_str();
  int na;
  get_data(masses_input,id_input,na,file_masses);
  
  vector<string> masses_req, id_req;
  
  // create WSTP link and compute amplitude
  if (!options.verbose)
  {
		create_wstp_link();
	}
  
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
    // need to get amplitude to find required masses
    print_math_body_1(input,options,get_cwd());
    send_to_math(input); // will return amplitude as a string
  }
  else
  {
    cout << "computing amplitude " << endl;
    templates::print_math_header(input);
    utils::print_math_body_1(input,options,get_cwd());
    
	  send_to_math(input); // will return amplitude as a string
	}
	
  const char* amplitude;
  if (!options.verbose)
  {
	  if(!WSGetString((WSLINK)pHandle, &amplitude))
	  {
	    cout << "Error getting amplitude from WSTP" << endl;
	  }
    
    std::pair <vector<string>,vector<string>> required = get_required_masses(masses_input,id_input, amplitude);
        
		masses_req = required.first;
		id_req = required.second;
        
    cout << "required masses are: ";
    
    for (unsigned int i = 0; i<masses_req.size(); i++)
    {
			cout << " " << masses_req[i];
		}
        
    cout << endl;
	}
	else
	{
		masses_req = masses_input;
		id_req = id_input;
	}
    
	if (check_if_available(options) && !options.force)
	{ 
	 // do nothing
	 utils::print_masses(input, masses_req);
	}
	else
	{ 
		utils::print_masses(input, masses_req);
    utils::print_math_body_2(input,options,masses_req);
    
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
  
  remove_fake_IR_divergence(input,masses_req);
  
  send_to_math(input);
  
  // split amplitude into order(1) and order(1/ma^2)
  
  if(std::find(masses_req.begin(), masses_req.end(), "ma") != masses_req.end()) 
  {
    input += "amp1 = Simplify[Coefficient[SelfEnergyFinite/. MassBuilderA[ma,4]-> (ma^2 MassBuilderA[ma,4]),ma,-2],TimeConstraint->100000];";
    seperate_amp = true;
	} else 
	{
    input += "amp1 = 0;";
	}
  
  input += "amp2 = Simplify[SelfEnergyFinite - amp1/ma^2,TimeConstraint->100000];";
  send_to_math(input);
  
  if (options.verbose)
  {
		return true;
	}
  
  // extract non-zero coefficients
  // send List of all possible basis integrals and get list back of coefficients
  
  // create std::map<std::string, Bases> containing all Bases objects
  full_basis = set_bases(masses_req, id_req);
  // obtain a corresponding vector<string> of unique Bases identifiers
  full_basis_id = extract_keys(full_basis);
  nb = full_basis_id.size();
  
  // generate a .m file to call within Mathematica (since it is a very long list)
  ofstream math_1;
  math_1.open (add_mpi_ext("output/math_1", options.mpi_process, "m"));
  print_math_basis(full_basis,math_1,"amp1","4","amp2");
  math_1 << "Export[\""<<get_cwd()<<"/output/output_"<< options.mpi_process << ".txt\", {" << endl;
  
  for (int i = 0; i < nb-1;i++)
  {
    math_1 << "{\""<<full_basis_id[i]<<" \", CForm[C1"<<full_basis_id[i]<<" /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
		math_1 << "{\""<<full_basis_id[i]<<" \", CForm[C2"<<full_basis_id[i]<<" /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
	}
  math_1 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C1"<<full_basis_id[nb-1]<<"  /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  math_1 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C2"<<full_basis_id[nb-1]<<"  /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
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
  for (int i=0; i<temp_int-1; i+=2)
  {
    full_basis[output_string[i]].coefficient1 = coeff_new[i];
    full_basis[output_string[i]].coefficient2 = coeff_new[i+1];
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
  print_math_basis(reduced_basis,math_2, "amp1","4","amp2");
  
  // construct the amplitudes as coefficients multiplied by basis integrals
  math_2 << "Amp1Trial = 0 ";
  for (int i = 0; i<nbr;i++)
  {
    math_2 << " + "<<reduced_basis_id[i]<< " * C1"<<reduced_basis_id[i]; // use the terms we know are non-zero
  }
  math_2 << ";"<<endl;
  
  math_2 << "Amp2Trial = 0 ";
  for (int i = 0; i<nbr;i++)
  {
    math_2 << " + "<<reduced_basis_id[i]<< " * C2"<<reduced_basis_id[i]; // use the terms we know are non-zero
  }
  math_2 << ";"<<endl;  
  
  // evaluate the difference between the trial and actual amplitude
  math_2 << "difference1 = Simplify[amp1-Amp1Trial];"<<endl;
  math_2 << "difference2 = Simplify[amp2-Amp2Trial];"<<endl;

  // evaluate the coefficients for the differences (see what basis integrals appear in the difference contains)
  // and write this to file
  print_math_basis(full_basis,math_2, "difference1",dimension,"difference2");
  math_2 << "Export[\""<<get_cwd()<<"/output/output2_"<< options.mpi_process << ".txt\", {" << endl;
  for (int i = 0; i < nb-1;i++)
  {
    math_2 << "{\""<<full_basis_id[i]<<" \", CForm[C1"<<full_basis_id[i]<<" + C2"<<full_basis_id[i]<<" + C1" <<full_basis_id[i] <<"2 + C2" <<full_basis_id[i] <<"2  /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
  }
  math_2 << "{\""<<full_basis_id[nb-1]<<" \", CForm[C1"<<full_basis_id[nb-1]<<" + C2"<<full_basis_id[nb-1]<<" + C1" <<full_basis_id[nb-1] <<"2 + C2" <<full_basis_id[nb-1] <<"2   /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
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
  
  print_math_basis(reduced_basis,math_3,"amp1","4","amp2");
  print_math_basis(prod_basis,math_3,"amp1","4","amp2");
  print_math_products(prod_basis,math_3,"amp1","4","amp2");
  
  math_3 << "SelfEnergyTrial1 = 0 ";
  // consider basis integrals that are not products
  for (int i = 0; i<nbr;i++)
  {
    math_3 << " + "<<reduced_basis_id[i]<< " * C1"<<reduced_basis_id[i];
  }
  // consider the products of basis integrals
  for (int i = 0; i<np ; i++)
  {
    for (int j = 0; j<np ; j++)
    {
      math_3 << " + "<< prod_id[i] << prod_id[j] << " * C1"<<prod_id[i]<< prod_id[j];
    }
  }
  math_3 << ";"<<endl;
  
  
   math_3 << "SelfEnergyTrial2 = 0 ";
  // consider basis integrals that are not products
  for (int i = 0; i<nbr;i++)
  {
    math_3 << " + "<<reduced_basis_id[i]<< " * C2"<<reduced_basis_id[i];
  }
  // consider the products of basis integrals
  for (int i = 0; i<np ; i++)
  {
    for (int j = 0; j<np ; j++)
    {
      math_3 << " + "<< prod_id[i] << prod_id[j] << " * C2"<<prod_id[i]<< prod_id[j];
    }
  }
  math_3 << ";"<<endl;
  
  
  // evaluate the difference between the constructed amplitude and the original amplitude
  math_3 << "diff1 = Simplify[amp1-SelfEnergyTrial1];\n";
  math_3 << "diff2 = Simplify[amp2-SelfEnergyTrial2];\n";

  // save the result to a text file
  math_3 << "Export[\""<<get_cwd()<<"/output/result1_"<< options.mpi_process << ".txt\", CForm[diff1/. DiracGamma[Momentum[p]] -> p] ]" << endl;
  math_3 << "Export[\""<<get_cwd()<<"/output/result2_"<< options.mpi_process << ".txt\", CForm[diff2/. Pair[Momentum[p], Momentum[p]]->p^2] ]" << endl;
  
  math_3 << "DumpSave[\""<<get_cwd()<<"/output/diff"<< options.mpi_process << ".mx\", diff2  ]" << endl;
  
  
  math_3 << "Export[\""<<get_cwd()<<"/output/output_products_"<< options.mpi_process << ".txt\", {" << endl;
  
  products_map.clear();
    
  // write out the final coefficients in C++ form
  for (int i = 0; i<np;i++)
  {
    for (int j = 0; j<np;j++)
    {
      math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C1"<<prod_id[i] << prod_id[j] <<"  /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
      if ((i==np-1) && (j==np-1))
      {
				 math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C2"<<prod_id[i] << prod_id[j] <<"  /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
			}
      else 
      {
				math_3 << "{\""<<prod_id[i] << prod_id[j] <<" \",CForm[C2"<<prod_id[i] << prod_id[j] <<" /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
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
  if (!options.verbose)
  {
		WSClose(link);
		cout << "WSTP link closed successfully" << endl;
	}
  
  
  
  ///// organise output data////////
  
  // check if the remainder contains any basis integrals
  string remainder1,remainder2;
  success = check_done(remainder1,remainder2,options.mpi_process);
  

  
  ReplaceAll(remainder1,"Pair(Momentum(p),Momentum(p))", "(-Power(p,2))");
  ReplaceAll(remainder2,"Pair(Momentum(p),Momentum(p))", "(-Power(p,2))");
  
  // BASIS INTEGRAL COEFFICIENTS //
  
  const char* coeff_integrals_tmp = "/output/coeff_integrals_";
  string c_coeff_integrals = "models/" + options.model +coeff_integrals_tmp + tag + ext;
  const char *coeff_integrals = c_coeff_integrals.c_str();
  ofstream coeff_integrals_out;
  coeff_integrals_out.open (coeff_integrals);
  format_coeff("4",reduced_basis,  reduced_basis_id, masses_req, id_req);
  if (remainder1 != "0")
  {
    format_coeff(remainder1,masses_req,id_req);
    coeff_integrals_out << "  TSIL_COMPLEXCPP C10 = " << remainder1 << ";" <<endl;
  }
  if (remainder2 != "0")
  {
    format_coeff(remainder2,masses_req,id_req);
    coeff_integrals_out << "  TSIL_COMPLEXCPP C20 = " << remainder2 << ";" <<endl;
  }  
  for (int i = 0; i < nbr;i++)
  {
    if (seperate_amp)
    {
			coeff_integrals_out << "  TSIL_COMPLEXCPP C1" << reduced_basis_id[i] << " = " << reduced_basis[reduced_basis_id[i]].coefficient1 << ";" <<endl;
		}
		coeff_integrals_out << "  TSIL_COMPLEXCPP C2" << reduced_basis_id[i] << " = " << reduced_basis[reduced_basis_id[i]].coefficient2 << ";" <<endl;
  }
  coeff_integrals_out.close();
  
  
  // PRODUCTS COEFFICIENTS //
  
  const char* file_integrals3_tmp = "output/output_products_";
  string c_file_integrals3 = file_integrals3_tmp + std::to_string(options.mpi_process) + ext;
  const char *file_integrals3 = c_file_integrals3.c_str();
  vector<string> name_products,coeff_products_new;
  
  get_data(name_products, coeff_products_new, temp_int,file_integrals3, true);
  std::map <std::string, Bases > prod_map = products_container(prod_id);
  
  for (int i=0; i<temp_int-1; i+=2)
  {
		prod_map[name_products[i]].coefficient1 = coeff_products_new[i];
    prod_map[name_products[i]].coefficient2 = coeff_products_new[i+1];
    
  }
  std::map <std::string, Bases > reduced_prod_map = remove_zeros(prod_map, extract_keys(prod_map));
  vector<string> reduced_prod_names = extract_keys(reduced_prod_map);
  
  const char* coeff_products_tmp = "/output/coeff_products_";
  string c_coeff_products = "models/" + options.model +coeff_products_tmp + tag + ext;
  const char *coeff_products = c_coeff_products.c_str();
  
  ofstream coeff_products_out;
  coeff_products_out.open (coeff_products);
  format_coeff("4",reduced_prod_map,  reduced_prod_names, masses_req, id_req);
  for (unsigned int i = 0; i < reduced_prod_names.size();i++)
  {
    if (seperate_amp)
    {
			coeff_products_out << "  TSIL_COMPLEXCPP C1" << reduced_prod_names[i] << " = " << reduced_prod_map[reduced_prod_names[i]].coefficient1 << ";" <<endl;
		}
		coeff_products_out << "  TSIL_COMPLEXCPP C2" << reduced_prod_names[i] << " = " << reduced_prod_map[reduced_prod_names[i]].coefficient2 << ";" <<endl;
  }
  coeff_products_out.close();
  

  //  SUMMATION
  
  
  const char* summation_tmp = "/output/summation_"; // vector containing file names
  string c_summation = "models/" + options.model +summation_tmp + tag + ext;
  const char *summation = c_summation.c_str();
  
  ofstream summation_out;
  summation_out.open (summation);
  
  if (seperate_amp)
  {
	  summation_out << "  TSIL_COMPLEXCPP result1 = ( ";
	  if (remainder1 != "0")
	  {
	    summation_out << " + C10 ";
	  }
	  
	  for (int i = 0; i<nbr;i++)
	  {
	    summation_out  << " + "<<reduced_basis_id[i]<< " * C1"<<reduced_basis_id[i];
	  }
	  for (unsigned int i = 0; i<reduced_prod_names.size();i++)
	  {
	    Bases temp_base;
	    temp_base = reduced_prod_map[reduced_prod_names[i]];
	    summation_out << " + "<< temp_base.e1 << " * " << temp_base.e2 << " * C1" << reduced_prod_names[i];
	  }
	  
	  if ((nbr==0) && (reduced_prod_names.size()==0) && (remainder1 == "0") && (remainder2 == "0"))
	  {
	    summation_out << "0";
	  }
	  
	  summation_out << ")/Power(ma,2);"<<endl;
	}
  
  
  summation_out << "  TSIL_COMPLEXCPP result2 = ";

  
  if (remainder2 != "0")
  {
    summation_out << " + C20 ";
  }
  
  
  for (int i = 0; i<nbr;i++)
  {
    summation_out  << " + "<<reduced_basis_id[i]<< " * C2"<<reduced_basis_id[i];
  }
  for (unsigned int i = 0; i<reduced_prod_names.size();i++)
  {
    Bases temp_base;
    temp_base = reduced_prod_map[reduced_prod_names[i]];
    summation_out << " + "<< temp_base.e1 << " * " << temp_base.e2 << " * C2" << reduced_prod_names[i];
  }
  
  if ((nbr==0) && (reduced_prod_names.size()==0) && (remainder1 == "0") && (remainder2 == "0"))
  {
    summation_out << "0";
  }
  
  summation_out << ";"<<endl;
  
  if (seperate_amp)
  {
		
		summation_out << "  if (exclude_photon_pole)\n"
		<< "  {\n"
		<< "    return result2;\n"
		<< "  }\n"
		<< "  else\n"
		<< "  {\n"
		<< "    return result1 + result2;\n"
		<< "  }\n";
	}
	else
	{
		summation_out << "  return result2;"<<endl;
	}
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
  
  if (options.verbose)
  {
		debug_out.open ("debug.m");
		debug_out << "(* ::Package:: *)" << endl;
		debug_out << "Quit[]" << endl;
		debug_out << "(* ::Section:: *)" << endl;
	}
	else
	{
		log_out.open ("output/log.m");
		log_out << "(* ::Package:: *)" << endl;
		log_out << "Quit[]" << endl;
		log_out << "(* ::Section:: *)" << endl;
	}
  
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
	string momentum = "Pair[Momentum[p],Momentum[p]]";
  if (options.fermion)
  {
		momentum = "p";
	}
  int nd = diagram.size();
  
  string dimension = "D";
  
  create_wstp_link();
  load_libraries();
  WSNewPacket(link);
  
  std::string input;
  

  input += "kappa=1/(16*Pi^2);";
  input += "SEtotal = 0 ;";
  
  for (int i=0;i<nd;i++)
  {
    
    input += "Get[\"" + get_cwd() + "/models/" + options.model + "/output/math_data_" + particle + "_" + diagram[i] + "_1" + ".mx\"];";
    
    input += "SEtotal = SEtotal + SelfEnergyFinite*kappa;";
  }
  
  
  input += "SE = makeFiniteAmplitude[SEtotal,-1,D];";
  
  // get the relevant counter-term amplitude
  
  input += "Get[\"" + get_cwd() + "/models/" + options.model + "/output/math_data_" + particle + "_1_1c" + ".mx\"];";
  
  input += "SEct = makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;";
  
  input += "ToString[SEct,CForm]";
  
  send_to_math(input);
  
  const char* counter_term_amplitude;
  
  if(!WSGetString((WSLINK)pHandle, &counter_term_amplitude))
  {
    cout << "Error getting string from WSTP" << endl;
  }
  
  // get the list of couplings from couplings.txt and find which ones are relevant here
  
  string c_file_couplings = "models/" + options.model + "/couplings.txt";  // need to make this model independent
  const char *file_couplings = c_file_couplings.c_str();
	
	int nc = 0;
	vector<std::string> couplings, relationships;
	vector<std::string> required_couplings;
  get_data(couplings, relationships, nc,file_couplings,true);

  string counter_term_amplitude_str = counter_term_amplitude;
  
  for (int i=0;i<nc;i++)
  {
	
		if (counter_term_amplitude_str.find(couplings[i]) != std::string::npos)
		{	
			stringstream _part_1;
			string part_1;
			_part_1 << couplings[i][0];
			_part_1 >> part_1;
      
			if (part_1=="d")
			{
				required_couplings.push_back(couplings[i]);
			}
		}
  }
	
	int nc_req = required_couplings.size();

	cout << "Required couplings are: ";
	for (int i = 0; i < nc_req; i++)
	{
		cout << " " << required_couplings[i];
	}
	cout << endl;
	
	// throw error if required couplings > 2

	// now solve relevant equations
	
	if (nc_req == 2)
	{
		input += "eq1 = FullSimplify[Coefficient[SE+SEct," + momentum + "]];";
		input += "eq2 = FullSimplify[Coefficient[SE+SEct," + momentum + ",0]];";
		input += "sol = Solve[{eq1==0,eq2==0},{" + required_couplings[0] + "," + required_couplings[1] + "}];";
		input += "Set @@@ sol[[1]];";
		send_to_math(input);
		
		// get each coupling
		
		input += "ToString[" + required_couplings[0] + ",CForm]";
		send_to_math(input);
  
		const char* coupling_1;
  
		if(!WSGetString((WSLINK)pHandle, &coupling_1))
		{
			cout << "Error getting string from WSTP" << endl;
		}
		cout << required_couplings[0] << " = " << coupling_1 << endl;
		input += "ToString[" + required_couplings[1] + ",CForm]";
		send_to_math(input);
  
		const char* coupling_2;
  
		if(!WSGetString((WSLINK)pHandle, &coupling_2))
		{
			cout << "Error getting string from WSTP" << endl;
		}
		cout << required_couplings[1] << " = " << coupling_2 << endl;		
		for (unsigned int i = 0; i < relationships.size() ; i++)
		{
			if (couplings[i] == required_couplings[0])
			{
				relationships[i] = coupling_1;
			}
			if (couplings[i] == required_couplings[1])
			{
				relationships[i] = coupling_2;
			}
		}
		
		
	}
	else if (nc_req == 1)
	{
		input += "eq1 = FullSimplify[SE+SEct];";
		input += "sol = Solve[{eq1==0},{" + required_couplings[0] + "}];";
		input += "Set @@@ sol[[1]];";
		input += "ToString[" + required_couplings[0] + ",CForm]";
		send_to_math(input);
  
		const char* coupling_1;
  
		if(!WSGetString((WSLINK)pHandle, &coupling_1))
		{
			cout << "Error getting string from WSTP" << endl;
		}
		cout << required_couplings[0] << " = " << coupling_1 << endl;
		for (unsigned int i = 0; i < relationships.size() ; i++)
		{
			if (couplings[i] == required_couplings[0])
			{
				relationships[i] = coupling_1;
			}
		}
		
		
	}
	else
	{
		cout << "could not determine 1 or 2 counter-term couplings to solve for." << endl;
	}
  
  
  input += "Quit[]";
  send_to_math(input);
  
  WSClose(link);
  cout << "WSTP link closed successfully" << endl;
  
  
  // now update the couplings list
  
  // we need couplings, relationships, required_couplings, coupling_1 and coupling_2
  
  ofstream output;
  output.open(file_couplings);
  
  
	for (unsigned int i = 0; i < relationships.size() ; i++)
	{
		output << couplings[i] << " " << trim_white_space(relationships[i]) << endl;
	}
	for (unsigned int i = relationships.size(); i < couplings.size(); i++)
	{
		output << couplings[i] << endl;
	}
}

void Compute_amp::calc_counter_terms()
{
  // need to read in the list of available diagrams and then select the 1-loop ones for
  // adding up here
  
  // read in available diagrams
  
  make_tag(options); // this is required to check for fermion
  open_log_files();
	
  
  const char *ext = ".txt";
  const char* file_diagrams_tmp = "models/";
  string c_file_diagrams = file_diagrams_tmp + options.model + "/output/avail_diagrams_" + ext;
  const char *file_diagrams = c_file_diagrams.c_str();
  
  if (options.input_list!="")
  {
		file_diagrams = options.input_list.c_str();
	}
  
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

