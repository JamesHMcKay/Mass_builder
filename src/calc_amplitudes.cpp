/*
Mass Builder - the missing link in automated two-loop self energy calculations
Please refer to the documentation for details or the readme.txt for simple run instructions

Last edited 29/09/16
James McKay

--- calc_amplitudes.cpp ---

generate mathematica scripts to compute all coefficients and find required basis integrals
also includes functions for sending FeynArts diagrams directory to pdf in the folder FAdiagrams

run ./mass_builder -g <model_name> after to generate code and ./mass_builder -e input.txt to evaluate

*/

#include "calc_amplitudes.hpp"

#define RUN_ALL
//#define DEBUG
using namespace std;
using namespace utils;

// create a mathematica object that is the product of two basis integrals

bool verbose=0;
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
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///   BEGIN NEW CONTENT ///
  
            vector<string> masses_input,id_input;
            int na;
            get_data(masses_input,id_input,na,file_masses);
            cout << "size of id input " << id_input.size() << endl;
            cout << "size of masses input " << masses_input.size() << endl;
            std::map<std::string, Bases> base_map = set_bases(masses_input, id_input);



            vector<string> bases_names = extract_keys(base_map);

            #ifdef DEBUG
            for (unsigned int i = 0; i < bases_names.size();i++)
            {
            Bases base_temp;
            base_temp = base_map[bases_names[i]];
            cout << bases_names[i] << endl;
            cout << "masses values are " << base_temp.e1 << " " << base_temp.e2 << " " << base_temp.e3 << " " << base_temp.e4 << " " << base_temp.e5 << " " << endl;
            }
            #endif

            // now need to print id = integral to myfile





            print_math_basis(base_map,myfile, "SEn");

            myfile << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;



            for (unsigned int i = 0; i < bases_names.size();i++)
            {
              myfile << "{\" "<<bases_names[i]<<" \", CForm[C"<<bases_names[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}," << endl;
            }

            myfile << "{\" "<<bases_names[bases_names.size()-1]<<" \", CForm[C"<<bases_names[bases_names.size()-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \"\"}" << endl;
            myfile << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
  
  
  
            myfile << "Export[\""<<s_cwd<<"/output/output2.txt\", {" << endl;



            for (unsigned int i = 0; i < bases_names.size();i++)
            {
              myfile << "{\"TSIL_COMPLEXCPP C"<<bases_names[i]<<" =\", CForm[C"<<bases_names[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
            }

            myfile << "{\"TSIL_COMPLEXCPP C"<<bases_names[bases_names.size()-1]<<" =\", CForm[C"<<bases_names[bases_names.size()-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;
            myfile << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;

  
  


  ///   END NEW CONTENT ///

  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  #ifdef RUN_ALL
  system("chmod +x output/stage_3.m ");
  if (verbose) system("./output/stage_3.m");
  else system("./output/stage_3.m  >/dev/null");
  system("chmod u+x scripts/script_1.sh");
  system("./scripts/script_1.sh");
  #endif
  
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////


  ///   BEGIN NEW CONTENT ///
  
  
  /*  new content below
  
  
  Idea here is to take output.txt and read in the coefficients to the Base object, so now each object in the map will have
  a coefficient associated with it.  We can then cull the map down to a much shorter map.  Then we will deal with products.
  
  In this section we reduce SEn trial to the cross terms
  
  */
  
  
  
  
            vector<std::string> output_string, coeff_new;
            int temp_int = 0;
            
            const char* file_integrals2_tmp = "output/output";
            string c_file_integrals2 = file_integrals2_tmp + blank + ext;
            const char *file_integrals2 = c_file_integrals2.c_str();
            
            get_data(output_string, coeff_new, temp_int,file_integrals2);
            for (int i=0; i<temp_int; i++)
            {
            base_map[output_string[i]].coefficient = coeff_new[i];
            }
            
            
            std::map <std::string, Bases > reduced_base_map = remove_zeros(base_map, bases_names);
            vector<string> reduced_bases_names = extract_keys(reduced_base_map);
            
            string prevb = "stage_3.mx\"";
            ofstream myfile_stage6b;
            myfile_stage6b.open ("output/stage_6b.m");
            
            utils::print_math_header(myfile_stage6b);
            myfile_stage6b<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n"
            <<"Print[SEn]"<< endl;
            

            
            print_math_basis(reduced_base_map,myfile_stage6b, "SEn");
            
            
            myfile_stage6b << "SEnTrial = 0 ";
            
            for (unsigned int i = 0; i<reduced_bases_names.size();i++)
            {
            
            myfile_stage6b << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
              
            }
            myfile_stage6b << ";"<<endl;
            
            myfile_stage6b << "diff = FullSimplify[SEn-SEnTrial]"<<endl;
            myfile_stage6b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage6b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage6b << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
            myfile_stage6b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage6b << "Print[\" --------------------------------------- \"]" <<endl;
            
            myfile_stage6b << "Print[diff]"<<endl;
            
            myfile_stage6b << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;
  
            print_math_basis(reduced_base_map,myfile_stage6b, "diff");
  
            myfile_stage6b << "Export[\""<<s_cwd<<"/output/output.txt\", {" << endl;


            if (reduced_bases_names.size()!=0)
            {
            for (unsigned int i=0;i<reduced_bases_names.size()-1;i++)
            {
            myfile_stage6b << "{\" TSIL_COMPLEXCPP C"<<reduced_bases_names[i]<<" =\", CForm[C"<<reduced_bases_names[i]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;
            }

            myfile_stage6b << "{\" TSIL_COMPLEXCPP C"<<reduced_bases_names[reduced_bases_names.size()-1]<<" =\", CForm[C"<<reduced_bases_names[reduced_bases_names.size()-1]<<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;
            myfile_stage6b << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;
            }

            myfile_stage6b << "Print[\"Completed\"]"<<endl;
  
            myfile_stage6b.close();

            
            
            
  ///  END NEW CONTENT ///
  
  
  
    //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  
  
  //////////////////////////////////////////////////////////////////
  //  Now create substatially more compact Mathematica routine
  //  using only the required basis integrals and non-zero coefficients
  //  then determine what cross-terms we still need to account for
  
  
  
  #ifdef RUN_ALL
      system("chmod +x output/stage_6b.m ");
      if(verbose) system("./output/stage_6b.m ");
      else system("./output/stage_6b.m  >/dev/null ");
  system("chmod u+x scripts/script_2.sh");
  system("./scripts/script_2.sh ");
  #endif
  
  
  n = 0;
  
//  if (check_done()) {success = 1;}
 // else
//  {
  
  
  //////////////////////////////////////////////////////////////////
  // of the terms that are remaining we now deal with possible combinations, the updated file for objects of interest
  // is names_updated.txt, this should contain all integrals for which we need to consider combinations
  
  
  ///  BEGIN NEW CONTENT ///
  
  
  

            vector<std::string> output_string_prod, coeff_new_prod;
            std::map <std::string, Bases > base_map_prod;
            vector<string> bases_names_prod;
            temp_int = 0;
  
            get_data(output_string_prod, coeff_new_prod, temp_int,file_integrals2); // file_integrals2 defined earlier, same file name
  
  
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
            reduced_base_map[output_string_prod[i]].coefficient = coeff_new_prod[i];
            }
            base_map_prod = remove_zeros(reduced_base_map, reduced_bases_names);
            bases_names_prod = extract_keys(base_map_prod);
            }

            ofstream myfile_stage8b;
            myfile_stage8b.open ("output/stage_8b.m");
            
            
            
            prevb = "stage_3.mx\"";
            utils::print_math_header(myfile_stage8b);
            myfile_stage8b<<"Get[\"" << s_cwd <<"/output/"<< prevb << "]\n"
            <<"Print[SEn]"<< endl;
  
            print_math_basis(reduced_base_map,myfile_stage8b,"SEn");
  
            print_math_products(base_map_prod,myfile_stage8b,"SEn");
  


            myfile_stage8b << "SEnTrial = ";
            for (unsigned int i = 0; i<reduced_bases_names.size();i++)
            {
            myfile_stage8b << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
            }




            for (unsigned int i = 0; i<bases_names_prod.size();i++)
            {
            for (unsigned int j = 0; j<bases_names_prod.size();j++)
            {
            myfile_stage8b << " + "<< bases_names_prod[i] << bases_names_prod[j] << " * C"<<bases_names_prod[i]<< bases_names_prod[j];
            }
            }
            myfile_stage8b << ";"<<endl;



  
            myfile_stage8b << "diff = FullSimplify[SEn-SEnTrial]"<<endl;
            myfile_stage8b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage8b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage8b << "Print[\" The difference between trial and actual SE is:\"]" <<endl;
            myfile_stage8b << "Print[\" --------------------------------------- \"]" <<endl;
            myfile_stage8b << "Print[\" --------------------------------------- \"]" <<endl;
            
            myfile_stage8b << "Print[diff]"<<endl;
            myfile_stage8b << "Export[\""<<s_cwd<<"/output/result.txt\", diff]" << endl;
  
            // print out coefficients of products
            int nn = bases_names_prod.size();
            myfile_stage8b << "Export[\""<<s_cwd<<"/output/output_products.txt\", {" << endl;

            for (int i = 0; i<nn;i++)
            {
            for (int j = 0; j<nn;j++)
            {
              
              if ((i==nn-1) && (j==nn-1)){ myfile_stage8b << "{\" TSIL_COMPLEXCPP C"<<bases_names_prod[i] << bases_names_prod[j] <<" =\", CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;}
              else {myfile_stage8b << "{\" TSIL_COMPLEXCPP C"<<bases_names_prod[i] << bases_names_prod[j] <<" =\", CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
              
            }
            }

            myfile_stage8b << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;

            // print out in a form useful for determing products later on

            myfile_stage8b << "Export[\""<<s_cwd<<"/output/output_products_2.txt\", {" << endl;

            for (int i = 0; i<nn;i++)
            {
            for (int j = 0; j<nn;j++)
            {
              
              if ((i==nn-1) && (j==nn-1)){ myfile_stage8b << "{\" "<<bases_names_prod[i] << "*" << bases_names_prod[j] << "*C"<<bases_names_prod[i] << bases_names_prod[j] <<"        =\", CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}" << endl;}
              else {myfile_stage8b << "{\" "<<bases_names_prod[i] << "*" << bases_names_prod[j] << "*C"<<bases_names_prod[i] << bases_names_prod[j] <<"           =\", CForm[C"<<bases_names_prod[i] << bases_names_prod[j] <<" /. Pair[Momentum[p], Momentum[p]] -> p^2 /. DiracGamma[Momentum[p]] -> p], \";\"}," << endl;}
              
            }
            }

            myfile_stage8b << " }, \"Table\", \"FieldSeparators\" -> \" \", \"TextDelimiters\" -> \"\"];" << endl;


            myfile_stage8b.close();




  
  
  
  
  ///  END NEW CONTENT ///
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  
  
  #ifdef RUN_ALL
  system("chmod +x output/stage_8b.m ");
  if(verbose) system("./output/stage_8b.m");
  else system("./output/stage_8b.m  >/dev/null ");
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

  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
  output_1 << reduced_bases_names[i] << endl;
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
  for (unsigned int i = 0; i<reduced_bases_names.size();i++)
  {
  if (sum_integrals != 0 ) summation_out  << " + "<<reduced_bases_names[i]<< " * C"<<reduced_bases_names[i];
  }
  
  const char* file_integrals_tmp = "output/names_updated";
  string c_file_integrals = file_integrals_tmp + blank + ext;
  const char *file_integrals = c_file_integrals.c_str();
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
