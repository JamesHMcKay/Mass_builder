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


#include "utils.hpp"

namespace utils
{
  time_t t = time(0);
  struct tm * now = localtime( & t );
  
  bool add_epsilon_terms = true;
  
  std::string get_cwd()
  {
    return getcwd(NULL,0);
  }
  
  
  void status_bar(double status)
  {
		std::string faces[6] = {":-C",":-(",":-|",":-)",":-D",":-O"};
    int fc = floor((status * 6.0)/100);
	  cout<< "\r" << "sorting integrals . . . " << status << "% complete " << faces[fc];
	  std::cout << std::flush;
	}
  
  
  
  int get_loop_order(string type)
  {
    int loop_order;
    string loop_order_str = utils::char_to_string(type[0]);
    stringstream convert(loop_order_str);
    convert >> loop_order;
    return loop_order;
  }
  
  // transform tag into suitable form
  // {a,b,c} -> group_a_b_c
  std::string clean(std::string input)
  {
		ReplaceAll(input,"{", "group_");
		ReplaceAll(input,",", "_");
		ReplaceAll(input,"}", "");
		return input;
	}
  
  
  const char * add_mpi_ext(std::string name, int process, std::string ext)
  {
    const char* file_tmp = "_";
    string c_file = name + file_tmp + std::to_string(process) + "." + ext;
    const char *file = c_file.c_str();
    return file;
  }
  
  
  const char * output_file_name(std::string model, std::string tag, std::string file)
  {
    const char *ext = ".txt";
    string underscore = "_";
    const char* str_tmp = "/output/";
    string c_str = "models/" + model +str_tmp + file + underscore + tag + ext;
    const char *str = c_str.c_str();
    return str;
  }
  
  
  void update_avail_diagrams(Options options)
  {
    const char *ext = ".txt";
    string underscore = "_";
    const char* str_tmp = "/output/";
    string c_str = "models/" + options.model +str_tmp + "avail_diagrams" + underscore + ext;
    const char *str = c_str.c_str();
    string file_name = str;
    
    vector<string> particle, diagram, type;
    vector<string> particle_new, diagram_new, type_new;
    int n;
    get_data(particle,diagram,type,n,file_name);
    
    ofstream file;
    file_name = output_file_name(options.model,"","avail_diagrams");
    file.open(file_name.c_str());
    
    particle.push_back(options.particle);
    diagram.push_back(options.diagram);
    
    if (options.counter_terms)
    {
      type.push_back( to_string(options.loop_order)+"c" ) ;
    }
    else
    {
      type.push_back( to_string(options.loop_order) ) ;
    }
    n = n+1;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if ((particle[j] == particle[i]) && (diagram[j]==diagram[i]) && (type[i]==type[j]) && (i!=j))
        {
          particle[i] = "";diagram[i]="";type[i]="";
        }
      }
    }
    int m = 0;
    for (int i = 0; i < n; i++)
    {
      if (particle[i]!="") { particle_new.push_back(particle[i]);diagram_new.push_back(diagram[i]);type_new.push_back(type[i]); m=m+1;}
    }
    for (int i = 0; i < m; i++)
    {
      file << particle_new[i] << " " << diagram_new[i] << " " << type_new[i] << endl;
    }
    file.close();
  }
  
  
  
  // check if a diagram has already been computed, return true if it has
  bool check_if_available(Options options)
  {
    const char *ext = ".txt";
    string underscore = "_";
    const char* str_tmp = "/output/";
    string c_str = "models/" + options.model +str_tmp + "avail_diagrams" + underscore + ext;
    const char *str = c_str.c_str();
    string file_name = str;
    
    vector<string> particle, diagram, type;
    int n;
    get_data(particle,diagram,type,n,file_name);
    
    string this_type;
    if (options.counter_terms)
    {
      this_type = to_string(options.loop_order)+"c";
    }
    else
    {
      this_type = to_string(options.loop_order);
    }
    
    for (int i = 0; i < n; i++)
    {
      if ((options.particle == particle[i]) && (options.diagram==diagram[i]) && (this_type==type[i]))
      {
        return true;
      }
    }
    
    return false;
  }
  
  
  
  
  void sort_avail_diagrams(Options options)
  {
    const char *ext = ".txt";
    string underscore = "_";
    const char* str_tmp = "/output/";
    string c_str = "models/" + options.model +str_tmp + "avail_diagrams" + underscore + ext;
    const char *str = c_str.c_str();
    string file_name = str;
    
    vector<string> particle, diagram, type;
    vector<string> particle_new, diagram_new, type_new;
    int n;
    get_data(particle,diagram,type,n,file_name);
    
    ofstream file;
    file_name = output_file_name(options.model,"","avail_diagrams");
    file.open(file_name.c_str());
    
    
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if ((particle[j] == particle[i]) && (diagram[j]==diagram[i]) && (type[i]==type[j]) && (i!=j))
        {
          particle[i] = "";diagram[i]="";type[i]="";
        }
      }
    }
    int m = 0;
    for (int i = 0; i < n; i++)
    {
      if (particle[i]!="") { particle_new.push_back(particle[i]);diagram_new.push_back(diagram[i]);type_new.push_back(type[i]); m=m+1;}
    }
    for (int i = 0; i < m; i++)
    {
      file << particle_new[i] << " " << diagram_new[i] << " " << type_new[i] << endl;
    }
    file.close();
  }
  
  
  void get_data(vector<std::string> &A,int &n,const char *filename)
  {
    n=0;
    std::ifstream input(filename);
    std::string line;
    while(getline(input, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      std::istringstream iss(line);
      n=n+1;
    }
    A.resize(n);
    input.close();
    n=0;
    std::ifstream input2(filename);
    std::string line2;
    while(getline(input2, line2))
    {
      if (!line2.length() || line2[0] == '#')
        continue;
      std::istringstream iss2(line2);
      iss2>> A[n];
      n=n+1;
    }
    input2.close();
  }
  
  
  
  std::string trim_white_space(std::string str)
	{
		size_t startpos = str.find_first_not_of(" \t");
		if( string::npos != startpos )
		{
			str = str.substr( startpos );
		}
		return str;
	}
  
  
  
  
  void get_data(vector<std::string> &A, vector<std::string> &B,int &n,const char *filename, bool whole_line)
  {
    n=0;
    std::ifstream input(filename);
    std::string line;
    int na=0,nb=0;
    while(getline(input, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      std::istringstream iss(line);
      std::string test;
      n=n+1;
      if (std::count( line.begin(), line.end(), ' ' )>0)
      {
        na=na+1;
        nb=nb+1;
      }
      else
      {
        na=na+1;
      }
    }
    
    A.resize(na);
    B.resize(nb);
    n = na;
    na=0,nb=0;
    std::ifstream input2(filename);
    std::string line2;
    while(getline(input2, line2))
    {
      if (!line2.length() || line2[0] == '#')
        continue;
      std::istringstream iss2(line2);
      if (std::count( line2.begin(), line2.end(), ' ' )>0)
      {
        if (whole_line)
        {
          iss2>> A[na]; getline(iss2,B[nb]);
          //B[nb] = trim_white_space(B[nb]);
        }
        else
        {
          iss2>> A[na] >> B[nb];
        }
        
        na=na+1;
        nb=nb+1;
      }
      else
      {
        iss2>> A[na];
        na=na+1;
      }
    }
    input.close();
    input2.close();
  }
  
  
  void get_data(vector<std::string> &A,vector<std::string> &B,vector<std::string> &C,int &n, string file_name_temp)
  {
    n=0;
    std::ifstream input(file_name_temp.c_str());
    std::string line;
    while(getline(input, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      std::istringstream iss(line);
      n=n+1;
    }
    input.close();
    A.resize(n);
    B.resize(n);
    C.resize(n);
    n=0;
    std::ifstream input2(file_name_temp.c_str());
    std::string line2;
    while(getline(input2, line2))
    {
      if (!line2.length() || line2[0] == '#')
        continue;
      std::istringstream iss2(line2);
      iss2>> A[n] >> B[n] >> C[n];
      n=n+1;
    }
    input2.close();
  }
  
  void get_data(vector<std::string> &A,vector<std::string> &B,vector<std::string> &C,vector<std::string> &D,int &n, string file_name_temp)
  {
    n=0;
    std::ifstream input(file_name_temp.c_str());
    std::string line;
    while(getline(input, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      std::istringstream iss(line);
      n=n+1;
    }
    input.close();
    A.resize(n);
    B.resize(n);
    C.resize(n);
    D.resize(n);
    n=0;
    std::ifstream input2(file_name_temp.c_str());
    std::string line2;
    while(getline(input2, line2))
    {
      if (!line2.length() || line2[0] == '#')
        continue;
      std::istringstream iss2(line2);
      iss2>> A[n] >> B[n] >> C[n] >> D[n];
      n=n+1;
    }
    input2.close();
  }
  
  
  
  // function to assign FCGV variables in default patched FeynArts models
  void assign_FCGV(ofstream &file,Options options)
  {
    // open file with conversion list
    const char* file_FCGV_tmp = "models/";
    string c_file_FCGV = file_FCGV_tmp + options.model + "/FCGV.txt";
    const char *file_FCGV = c_file_FCGV.c_str();
    
    vector<string> variable,replacement;
    int n;
    get_data(variable,replacement,n,file_FCGV);
    
    for (int i = 0; i < n ;i++)
    {
      file << "FCGV[\"" <<variable[i] << "\"] = " << replacement[i] << ";";
    }
  }
  
  void assign_FCGV(std::string &input,Options options)
  {
    // open file with conversion list
    const char* file_FCGV_tmp = "models/";
    string c_file_FCGV = file_FCGV_tmp + options.model + "/FCGV.txt";
    const char *file_FCGV = c_file_FCGV.c_str();
    
    vector<string> variable,replacement;
    int n;
    get_data(variable,replacement,n,file_FCGV);
    
    for (int i = 0; i < n ;i++)
    {
      input+= "FCGV[\"" + variable[i] + "\"] = " + replacement[i] + ";";
    }
  }
  
  // function to reassign variable names, such as mixing matrices
  void assign_variables(ofstream &file,Options options)
  {
    // open file with conversion list
    const char* file_var_tmp = "models/";
    string c_file_var = file_var_tmp + options.model + "/reassign_variables.txt";
    const char *file_var = c_file_var.c_str();
    
    vector<string> variable,replacement;
    int n;
    get_data(variable,replacement,n,file_var);
    
    for (int i = 0; i < n ;i++)
    {
      file <<  variable[i] << " = " << replacement[i] << ";";
    }
    
  }
  
  void assign_variables(std::string &input,Options options)
  {
    // open file with conversion list
    const char* file_var_tmp = "models/";
    string c_file_var = file_var_tmp + options.model + "/reassign_variables.txt";
    const char *file_var = c_file_var.c_str();
    
    vector<string> variable,replacement;
    int n;
    get_data(variable,replacement,n,file_var);
    
    for (int i = 0; i < n ;i++)
    {
      input+=  variable[i] + " = " + replacement[i] + ";";
    }
    
  }
  
  string make_tag(Options &options)
  {
    // return safe name for input and output particles
    string particle_1 =  part_name_simple(options.particle_1);
    string particle_2 =  part_name_simple(options.particle_2);
    
    // define tag used in all generated files
    string particle_tag;
    if ((options.particle_1!=options.particle_2))
    {
      particle_tag = particle_1 + "_" + particle_2;
    }
    else
    {
      particle_tag = particle_1;
    }
    string tag="";
    
    if (options.counter_terms)
    {
      tag = particle_tag + "_" + options.diagram + "_" + to_string(options.loop_order)+"c";
    }
    else
    {
      tag = particle_tag + "_" + options.diagram + "_" + to_string(options.loop_order);
    }
    
    
    stringstream _part_1;
    string part_1;
    _part_1 << options.particle[0];
    _part_1 >> part_1;
      
    if (part_1=="F")
    {
      options.fermion = true;
    }
    if (part_1=="V")
    {
      options.vector = true;
    }
    
    return tag;
  }
  
  
  void get_finite_amplitude(std::string &input, Options options)
  {
    
    // obtain finite peice of amplitude (using MassBuilder package)
    if (options.counter_terms)
    {
      input += "SelfEnergyFinite = makeFiniteCT[SelfEnergyFinite," + options.epsilon_order + ", D];";
    }
    else
    {
      input += "SelfEnergyFinite = makeFiniteAmplitude[SelfEnergyFinite," + options.epsilon_order + ", D];";
    }
    
    // take transverse part of self energy (fix required after adding Truncated->True to get diagrams, as this removes
    // the polarization vectors from the amplitude which would normally fix this problem (now only used for vectors)
    input += "SelfEnergyFinite = SelfEnergyFinite/.(Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] -Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> Pair[Momentum[p],Momentum[p]] ;";
    input += "SelfEnergyFinite = SelfEnergyFinite/.(-Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] +Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> -Pair[Momentum[p],Momentum[p]] ;";
    
    input += "SelfEnergyFinite = Simplify[SelfEnergyFinite/.MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]] /. MassBuilderP -> p /. MassBuilderQ2->Q2 /. MassBuilderZeta-> Zeta,TimeConstraint->100000];";

  }
  
 
  void get_saved_amplitude(std::string &input, Options options, std::string mass)
  {
    
    //assign_FCGV(input,options);
    
    //assign_variables(input,options);
   
    string tag = make_tag(options);
    
    // define kinematics
    if (options.onshell && (options.loop_order!=1))
    {
      input+="SPD[p, p] = " + mass + "^2;";
      input+="Pair[Momentum[p], Momentum[p]] = " + mass + "^2;";
    }
    
    input += "Get[\"" + get_cwd() + "/models/" + options.model + "/output/math_data_" + tag + ".mx\"];";
  }
  
  
  
  
  void print_math_body_1(std::string &input,Options options,string cwd)
  {
    int loop_order = options.loop_order;
    string particle_1 = options.particle_1;
    string particle_2 = options.particle_2;
    string diagram = options.diagram;
    string model = options.model;
    
    assign_FCGV(input,options);
    assign_variables(input,options);
    
    
    
    if (options.counter_terms)
    {
      input+="t12 = CreateCTTopologies[" + to_string(loop_order) + ", 1 ->  "  +  to_string(options.n_final_states) + ", ExcludeTopologies -> Internal];";
    }
    else
    {
      input+="t12 = CreateTopologies[" +  to_string(loop_order) + ", 1 -> "  +  to_string(options.n_final_states)  +  ", ExcludeTopologies -> Internal];";
    }
    
    if (options.use_lorentz)
    {
      input+="alldiags = InsertFields[t12, {" + particle_1 + "} -> {" + particle_2 + "},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions -> {"  +  options.restrictions  +  "},Model -> \"" + cwd + "/models/" + model + "/" + model + "\"];";
    }
    else
    {
      input+="alldiags = InsertFields[t12, {" + particle_1 + "} -> {" + particle_2 + "},InsertionLevel -> {Particles}, GenericModel -> \"" + cwd + "/models/" + model + "/" + model + "\",Restrictions -> {"  +  options.restrictions  +  "},Model -> \"" + cwd + "/models/" + model + "/" + model + "\"];";
    }
    
    input+="subdiags0 =   DiagramExtract[alldiags, " + diagram + "];";
    
    
    string truncated = "True";
    if (options.vector){ truncated = "False";}
    
    input+="amp0 = FCFAConvert[";
    input+="CreateFeynAmp[subdiags0,Truncated -> " + truncated;
    
    if (options.general_gauge)
    {
			input+=", GaugeRules -> {}";
		}
    input+="], ";
    input+="IncomingMomenta -> {p}, OutgoingMomenta -> {p}, LoopMomenta -> {k1, k2} ,";
    input+="UndoChiralSplittings -> True,TransversePolarizationVectors -> {p}, ";
    input+="DropSumOver -> False, List -> False,ChangeDimension -> D] // Contract // FCTraceFactor;";
    input+="amp0 = amp0 /.SumOver[SUNFIndex[x_], 3]->3;";
    if (options.general_gauge)
    {
			input+="amp0 = amp0 /. GaugeXi[Z] -> Xi^2 /. GaugeXi[A] -> Xi^2 /. GaugeXi[W] -> Xi^2 /. GaugeXi[P] -> Xi^2 /. GaugeXi[Wp] -> Xi^2; ";
    }
    
    // correct for additional factors of Pi introduced by FeynArts
    
    if (options.counter_terms)
    {
      if (options.loop_order == 2 )
      {
        input+= "amp0 = amp0*16*Pi^4;";
      }
    }
    else
    {
      if (options.loop_order == 1 )
      {
        input+= "amp0 = amp0*16*Pi^4;";
      }
      if (options.loop_order == 2 )
      {
        input+= "amp0 = amp0*256*Pi^8;";
      }
    }
  
		input += "ToString[amp0/.Sqrt[Xi^2]->Xi,InputForm]";
  
  }


	void print_masses(std::string &input, std::vector<std::string> masses)
	{
		input+= "masses = List["  +  masses[0];
    if (masses.size()>1)
    {
      for (unsigned int i=1; i < ( masses.size() ) ;i++)
      {
        input+= ","  +  masses[i];
      }
    }
    input+= "];";
    
    if (masses[0]!="null")
    {
      input+= "massesExpand = List[1";
    }
    else
    {
      input+= "massesExpand = List[0";
    }
  
    if (masses.size()>1)
    {
      for (unsigned int i=1; i < ( masses.size() ) ;i++)
      {
        if (masses[i]!="null")
        {
          input+= ",1";
        }
        else
        {
          input+= ",0";
        }
      }
    }
    input+= "];\n";
  }




	void print_math_body_2(std::string &input,Options options,std::vector<std::string> masses)    
	{
		int loop_order = options.loop_order;
    string particle_1 = options.particle_1;
    string particle_2 = options.particle_2;
    string diagram = options.diagram;
    string model = options.model;
    
    
    
    for (int index = 1; index < 4;index ++)
    {
      input+="amp0 = amp0 /. Index[Generation, " + to_string(index) + "] -> " + to_string(index)+ ";";
    }

        
    if ( (loop_order == 2) && (!options.counter_terms) )
    {
      input += "ampsSE1 = (amp0 /. DiracTrace -> Tr) // FCMultiLoopTID[#, {k1, k2}] & // DiracSimplify;";
    }
    else if ( (loop_order == 1) && (options.counter_terms) )
    {
      input += "ampsSE1 = amp0;";
    }
    else
    {
      input+=" ampsSE1 = (amp0) // DiracSimplify // TID[#, k1] & // DiracSimplify;";
    }
    
    // define kinematics
    if (options.onshell && (options.loop_order!=1))
    {
      input+="SPD[p, p] = " + masses[0] + "^2;";
      input+="Pair[Momentum[p], Momentum[p]] = " + masses[0] + "^2;";
    }
    
    
    if (options.fire)
    {
      if ( (loop_order == 2) && (!options.counter_terms) )
      {
        input += "ampsSE2 = ampsSE1 // FDS[#, k1, k2] &;";
        input += "ampsSE3 = FIREBurn[ampsSE2, {k1, k2}, {p}] // FDS[#, k1, k2] &;";
        input += "ampsSE4 = ampsSE3 // Collect2[#, {FeynAmpDenominator}, Factoring -> Simplify] &;";
        input += "resSE = Cancel[ampsSE4];";
        input += "tfiamp0 = resSE // ToTFI[#, k1, k2, p] &;";
      }
      
      else if ( (loop_order == 1) && (options.counter_terms) )
      {
        input+=" fullamp0 = (ampsSE1) // DiracSimplify;";
        input+="tfiamp0 = fullamp0 // ChangeDimension[#, D] &;";
      }
      else
      {
        input += "ampsSE3 = FIREBurn[ampsSE1, {k1}, {p}] // FDS[#, k1] &;";
        input += "ampsSE4 = ampsSE3 // Collect2[#, {FeynAmpDenominator}, Factoring -> Simplify] &;";
        input += "resSE = Cancel[ampsSE4];";
        input += "tfiamp0 = resSE // ToTFI[#, k1, p] &;";
      }
    }
    else
    {
      input+="amp1 = ampsSE1 //FDS[#,k1,k2]&;";
      
      input+="SetOptions[Eps, Dimension -> D];";
      
      if ( (loop_order == 2) && (!options.counter_terms) )
      {
        input+="fullamp0 = (amp1) // DiracSimplify // FCMultiLoopTID[#, {k1, k2}] & //DiracSimplify;";
        input+="tfiamp0 = fullamp0 // ToTFI[#, k1, k2, p] & // ChangeDimension[#, D] &;";
      }
      
      else if ( (loop_order == 1) && (options.counter_terms) )
      {
        input+=" fullamp0 = (amp1) // DiracSimplify;";
        input+="tfiamp0 = fullamp0 // ChangeDimension[#, D] &;";
      }
      else
      {
        input+=" fullamp0 = (amp1) // DiracSimplify // TID[#, k1] & // DiracSimplify;";
        input+="tfiamp0 = fullamp0 // ToTFI[#, k1, p] & // ChangeDimension[#, D] &;";
      }
    }



  }
  
  void print_tarcer_recurse(std::string &input,Options options)
  {
    input+="SE = TarcerRecurse[tfiamp0];";
    input+="SEk = (1/(4 Pair[Momentum[p, D],Momentum[p, D]])) DiracTrace[ DiracGamma[Momentum[p, D], D].SE ];";
    input+="SEm = (1/4) DiracTrace[ SE ];";
    input+="SE = p*SEk+SEm;"; // TODO make this general (this only works for vectors because SEk is just zero and the p^2 comes through SEm)
    input+="SE = SE /. Pair[Momentum[Polarization[p, -I, Transversality -> True], D], Momentum[Polarization[p, I, Transversality -> True], D]] -> -1 ;";
    input+="SE = SE /. Sqrt[Xi^2]->Xi;"; // avoid doing a simplify (with specified Xi > 0) just to change this
    // uncomment the following line if using a different gauge choice
    if (options.general_gauge)
    {
			input+="SEn = SEn /. GaugeXi[Z] -> 0 /. GaugeXi[P] -> 0 /. GaugeXi[Wp] -> 0  /. GaugeXi[S[1]] -> 0 /. GaugeXi[S[2]] -> 0 /. GaugeXi[S[3]] -> 0 ;\n";
		}
  }
  
  
  void remove_fake_IR_divergence(std::string &input, std::vector<std::string> masses)
  {
		std::string massesSmall[10] = {"ma","mf","md","mu","ms","mb","mc","mm","ml","me"};
		bool small_mass = false;
    
		input+= "massesSmall = List[";
		
		for (unsigned int i = 0; i < masses.size(); i++)
		{
			small_mass = false;
			for (int j = 0; j < 10 ; j++)
			{
				if (masses[i] == massesSmall[j])
				{
					small_mass = true;
				}
			}
			
			if ( small_mass && (i < masses.size()-1))
			{
				input+= "1,";
			}
			else if (i < masses.size()-1)
			{
				input+= "0,";
			}
		}
		
		if (small_mass)
		{
			input+= "1];";
		}
		else
		{
			input+= "0];";
		}
		
		input+= "SelfEnergyFinite = Simplify[implementTbar[SelfEnergyFinite,";
		input+= "masses, massesSmall,MassBuilderA, MassBuilderB, MassBuilderT],TimeConstraint->100000];";
		
		// do a test with and without this for F11 34
		//input+= "SelfEnergyFinite = Limit[SelfEnergyFinite, ma->0, Analytic -> True ];";
	}
  
  
  
  bool ContainsXi(std::string input)
  {
    bool containsXi = false;
		if (input.find("Xi") != std::string::npos)
		{
			containsXi = true;
		}
    return containsXi;
  }
  
  
  std::vector<std::string> remove_xi(std::vector<std::string> input)
  {
		int nm = input.size();
		std::string from = "*";
		std::string to = "";
		std::vector<string> output;
		for (int i = 0; i < nm ; i++)
		{
			if (!ContainsXi(input[i]))
			{
				output.push_back(input[i]);
			}
		}
		
		return output;
	}
  
  
  std::vector<std::string> remove_multiplication(std::vector<std::string> input)
  {
		int nm = input.size();
		std::string from = "*";
		std::string to = "";
		for (int i = 0; i < nm ; i++)
		{
			ReplaceAll(input[i],from, to);
		}
		
		return input;
	}
  
  // find the masses that appear in the amplitude (given a list of options)
  std::pair <vector<string>,vector<string>> get_required_masses(vector<string> masses_in, vector<string> id_in, string amp)
  {
		int nm = masses_in.size();
		vector<string> masses_out;
		vector<string> id_out;
		
		for (int i = 0; i<nm; i++)
    {
			string str = " " + masses_in[i] + "]";
			if (masses_in[i] == "null" )
			{
				str = " 0]";
			}
			
      if (amp.find(str) != std::string::npos)
      {
        masses_out.push_back(masses_in[i]);
        id_out.push_back(id_in[i]);
      }
    }
    
    if (masses_out.size()==0)
    {
			masses_out.push_back("none");
			id_out.push_back("n");
		}
    
    std::pair <vector<string>,vector<string>> result;
    
    result.first = masses_out;
    result.second = id_out;
    
    return result;
    
	}
  
  bool check_done_quiet(int mpi_process)
  {
    std::ifstream file1(add_mpi_ext("output/result1",mpi_process,"txt"));
    std::ifstream file2(add_mpi_ext("output/result2",mpi_process,"txt"));
    std::string str1,str2;
    std::string result;
    std::getline(file1, str1);
    std::getline(file2, str2);
    result += str1 + str2;
    // need to check if the result contains a basis integral
    bool success = true;
    
    std::string bad_list[6] = {"TFI","TBI","TVI","TAI","TJI","MassBuilderTBAR"};
    
    for (int i = 0; i<6; i++)
    {
      if (result.find(bad_list[i]) != std::string::npos)
      {
        success = false;
      }
    }
    
    
    return success;
  }
  
  bool check_done(string &result1,string &result2, int mpi_process)
  {
		ifstream file1(add_mpi_ext("output/result1",mpi_process,"txt"));
    ifstream file2(add_mpi_ext("output/result2",mpi_process,"txt"));
    string str1,str2; 
    getline(file1, str1);
    getline(file2, str2);
    string result = str1 + str2;
    result1 += str1;
    result2 += str2;
    // need to check if the result contains a basis integral
    bool success = true;
    
    std::string bad_list[6] = {"TFI","TBI","TVI","TAI","TJI","MassBuilderTBAR"};
    
    for (int i = 0; i<6; i++)
    {
      if (result.find(bad_list[i]) != std::string::npos)
      {
        success = false;
      }
    }
    
    
    if (success)
    {
      cout << "Successful!!!" << endl;
    }
    else
    {
      cout << "Something has gone wrong.  Check that all mass terms match the masses in the model file and that there are no masses missing in your masses.txt input file." << endl;
      cout << "Also confirm that Mathematica is starting correctly and finding the required packages, run with \"-v\" flag to display Mathematica output." << endl;
    }
    return success;
  }
  
  // TODO: make this routine more generic for any kind of particle name and number
  std::string part_name_simple(std::string particle_name_full)
  {
    stringstream _part_1,_part_2;
    string part_1,part_2;
    
    if (particle_name_full.size() == 4)
    {
      _part_1 << particle_name_full[0];
      _part_1 >> part_1;
      _part_2 << particle_name_full[2];
      _part_2 >> part_2;
      return part_1+part_2;
    }
    else if (particle_name_full.size() == 9)
    {
      stringstream _part_3,_part_4;
      string part_3,part_4;
      
      _part_1 << particle_name_full[0];
      _part_1 >> part_1;
      _part_2 << particle_name_full[2];
      _part_2 >> part_2;
      _part_3 << particle_name_full[3];
      _part_3 >> part_3;
      
      _part_4 << particle_name_full[6];
      _part_4 >> part_4;
      
      return part_1+part_2+part_3+"_g"+part_4;
    }
    else
    {
      cout << "An error has occured in determing a safe name for the particles, only particle names X[n] or X[nn,{n}] are currently supported, please add required case" << endl;
      return particle_name_full;
    }
    
  }
  
  std::string part_name_simple(std::string particle_name_full_1,std::string particle_name_full_2)
  {
    std::string simple_1 = part_name_simple(particle_name_full_1);
    std::string simple_2 = part_name_simple(particle_name_full_2);
    
    if (simple_1==simple_2)
    {
      return simple_1;
    }
    else
    {
      return simple_1 + "_" + simple_2;
    }
  }
  
  
  vector<string> remove_duplicates(vector<string> input,string warning)
  {
    vector<std::string> input_unique = input;
    sort(input_unique.begin(),input_unique.end());
    input_unique.erase( unique( input_unique.begin(), input_unique.end() ), input_unique.end() );
    if (input.size() != input_unique.size())
    {
      cout << warning <<endl;
    }
    return input_unique;
  }
  
  vector<string> remove_duplicates(vector<string> input)
  {
    vector<std::string> input_unique = input;
    sort(input_unique.begin(),input_unique.end());
    input_unique.erase( unique( input_unique.begin(), input_unique.end() ), input_unique.end() );
    return input_unique;
  }
  
  vector<char> remove_duplicates(vector<char> input,string warning)
  {
    vector<char> input_unique = input;
    sort(input_unique.begin(),input_unique.end());
    input_unique.erase( unique( input_unique.begin(), input_unique.end() ), input_unique.end() );
    if (input.size() != input_unique.size())
    {
      cout << warning <<endl;
    }
    return input_unique;
  }
  
  vector<char> remove_duplicates(vector<char> input)
  {
    vector<char> input_unique = input;
    sort(input_unique.begin(),input_unique.end());
    input_unique.erase( unique( input_unique.begin(), input_unique.end() ), input_unique.end() );
    return input_unique;
  }
  
  
  
  string char_to_string(char c)
  {
    stringstream ss;
    string s;
    ss << c;
    ss >> s;
    return s;
  }
  
  
  vector<int> find_string_lengths(vector<string> input)
  {
    vector<int> lengths;
    for (unsigned int i = 0; i < input.size(); i++)
    {
      string input_temp = input[i];
      lengths.push_back( input_temp.size() );
    }
    return lengths;
  }
  
  
  std::vector<std::string> extract_keys(std::map<std::string, Bases> const& input_map)
  {
    std::vector<std::string> retval;
    for (auto const& element : input_map)
    {
      retval.push_back(element.first);
    }
    return retval;
  }
  
  
  void print_base(ofstream &myfile, Bases base, string id, string target, string D, string target2)
  {
    string type = base.type;
    string momentum = "Pair[Momentum[p],Momentum[p]]";
    if (D=="D")
    {
      momentum = "Pair[Momentum[p,D],Momentum[p,D]]";
    }
    
    if (target != "diff")
    {
      if (type == "F")
      {
        myfile << id << " = " << "TFI["<<D<<", " << momentum << ", {{1, " << base.e1  << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}, {1, " << base.e4 << "}, {1, " << base.e5 << "}}];" << endl;
      }
      if (type == "A" && base.e1!="null")
      {
        myfile << id << " = " << "TAI["<<D<<", 0, {{1, " << base.e1 << "}}];" << endl;
      }
      if (type == "B")
      {
        myfile << id << " = " << "TBI["<<D<<", " << momentum << ", {{1, " << base.e1 << "}, {1, " << base.e2 << "}}];" << endl;
      }
      if (type == "V")
      {
        myfile << id << " = " << "TVI["<<D<<", " << momentum << ", {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}, {1, " << base.e4 << "}}];" << endl;
      }
      if (type == "T")
      {
        myfile << id << " = " << "TJI["<<D<<", " << momentum << ", {{2, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      }
      if (type == "TBAR")
      {
        myfile << id << " = " << "MassBuilderTBAR[" << base.e1 << ", " << base.e2 << ", " << base.e3 << "];" << endl;
      }
      if (type == "J")
      {
        myfile << id << " = " << "TJI["<<D<<", " << momentum << ", {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      }
      if (type == "K" && (base.e1!="null" || base.e2!="null" || base.e3!="null"))
      {
        myfile << id << " = " << "TJI["<<D<<", 0, {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      }
    }
    
    myfile << "C1"<< id << " = Simplify[Coefficient["<<target<<", " << id << ", 1],TimeConstraint->100000];" << endl;
    myfile << "C1"<< id << "2 = Simplify[Coefficient["<<target<<", " << id << ", 2],TimeConstraint->100000];" << endl;
    
    if (target2!="")
    {
			myfile << "C2"<< id << " = Simplify[Coefficient["<<target2<<", " << id << ", 1],TimeConstraint->100000];" << endl;
			myfile << "C2"<< id << "2 = Simplify[Coefficient["<<target2<<", " << id << ", 2],TimeConstraint->100000];" << endl;
	  }
    
  }
  
  // print the basis integrals out in Mathematica notation
  void print_math_basis(std::map<std::string, Bases> base_map, ofstream &myfile, string target, string D,string target2)
  {
    vector<string> bases_names = extract_keys(base_map);
    for (unsigned int i = 0; i < bases_names.size();i++)
    {
      Bases base_temp;
      base_temp = base_map[bases_names[i]];
      print_base(myfile, base_temp, bases_names[i], target, D, target2);
    }
  }
  
  void print_base_product(ofstream &myfile,Bases base_1,Bases base_2,string target, string D, string target2)
  {
    string type1 = base_1.type;
    string type2 = base_2.type;
    string momentum = "Pair[Momentum[p],Momentum[p]]";
    if (D=="D")
    {
      momentum = "Pair[Momentum[p,D],Momentum[p,D]]";
    }
    if (type1 =="F"){goto end;}
    if (type2 =="F"){goto end;}
    myfile << base_1.short_name << base_2.short_name << " = ";
    if (type1=="A") myfile << "TAI["<<D<<", 0, {{1, " << base_1.e1 << "}}]";
    if (type1=="B") myfile << "TBI["<<D<<", " << momentum << " , {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}}]";
    if (type1=="J") myfile << "TJI["<<D<<", " << momentum << ", {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}}]";
    if (type1=="T") myfile << "TJI["<<D<<", " << momentum << ", {{2, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}}]";
    if (type1=="TBAR") myfile << "MassBuilderTBAR[" << base_1.e1 << ", " << base_1.e2 << ", " << base_1.e3 << "]";
    if (type1=="K") myfile << "TJI["<<D<<", 0, {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}}]";
    if (type1=="V") myfile << "TVI["<<D<<", " << momentum << ", {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "}}]";
    //if (type1=="F") myfile << "TFI["<<D<<", Pair[Momentum[p,D],Momentum[p,D]], {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "},{1, " << base_1.e5 << "}}]";
    myfile << " * ";
    if (type2=="A") myfile << "TAI["<<D<<", 0, {{1, " << base_2.e1 << "}}];";
    if (type2=="B") myfile << "TBI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}}];";
    if (type2=="J") myfile << "TJI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];";
    if (type2=="T") myfile << "TJI["<<D<<", " << momentum << ", {{2, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];";
    if (type2=="TBAR") myfile << "MassBuilderTBAR[" << base_2.e1 << ", " << base_2.e2 << ", " << base_2.e3 << "];";
    if (type2=="K") myfile << "TJI["<<D<<", 0, {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];";
    if (type2=="V") myfile << "TVI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}, {1, " << base_2.e4 << "}}];";
    //if (type2=="F") myfile << "TFI["<<D<<", Pair[Momentum[p,D],Momentum[p,D]], {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "},{1, " << base_1.e5 << "}}]";
    myfile << endl;
    
    
    if (base_1.short_name==base_2.short_name) myfile << "C1"<< base_1.short_name << base_2.short_name << " = Coefficient["<<target<<","<< base_1.short_name << ", 2];" << endl;
    else myfile << "C1"<< base_1.short_name << base_2.short_name << " = - (1/2)* Coefficient["<<target<<","<< base_1.short_name << base_2.short_name << ", 1];" << endl;
		if (target2!="")
		{
			if (base_1.short_name==base_2.short_name) myfile << "C2"<< base_1.short_name << base_2.short_name << " = Coefficient["<<target2<<","<< base_1.short_name << ", 2];" << endl;
			else myfile << "C2"<< base_1.short_name << base_2.short_name << " = - (1/2)* Coefficient["<<target2<<","<< base_1.short_name << base_2.short_name << ", 1];" << endl;
		}
  
  end:;
  }
  
  
  void print_math_products(std::map<std::string, Bases> base_map, ofstream &myfile, string target, string D, string target2)
  {
    vector<string> bases_names = extract_keys(base_map);
    int n = bases_names.size();
    for (int i = 0; i<n;i++)
    {
      for (int j = 0; j<n;j++)
      {
        
        Bases base_1;
        base_1 = base_map[bases_names[i]];
        base_1.short_name = bases_names[i];
        Bases base_2;
        base_2 = base_map[bases_names[j]];
        base_2.short_name = bases_names[j];  // TODO make this an automatic get and set function perhaps, could do short_names nicer
        
        print_base_product(myfile,base_1,base_2,target,D,target2);
        
      }
    }
  }
  
  void ReplaceAll(std::string &input, const std::string& from, const std::string& to)
  {
    size_t start_pos = 0;
    while((start_pos = input.find(from, start_pos)) != std::string::npos)
    {
      input.replace(start_pos, from.length(), to);
      start_pos += to.length();
    }
  }
  
  void print_doTSIL(ofstream &myfile,Bases base)
  {
    string type = base.type;
    string name = base.short_name;
    
    if (type == "const" || type == "A" || type == "B")
    {
			// do nothing
		}
		else
		{
			myfile << " if (data.do_tsil_all) " << endl;
			myfile << " { " << endl;
		}
    
    if (type == "const")
    {
      myfile << "  const" << " = 1.0L"<<endl;
    }
    if (type == "A")
    {
      myfile << "  " <<  name << " = -i*TSIL_A_ ("<<base.e1<<"2 , Q2);"<<endl;
    }
    
    set_ma_zero(base);
    
    if (type == "B")
    {
      myfile << "  " << name <<" = i*TSIL_B_ (" << base.e1 << "2, " << base.e2 << "2, s, Q2);"<< endl;
    }
    if (type == "K")
    {
      myfile << "  " << name <<" = TSIL_I2_(" << base.e1 << "2, " << base.e2 << "2, " << base.e3 << "2, Q2);"<< endl;
    }
    if (type == "J")
    {
      string A = base.e1,B=base.e2,C=base.e3;
      myfile << "  TSIL_SetParametersST (&bar," << B << "2, " << A << "2, " << C <<"2, Q2);" << endl;
      myfile << "  TSIL_Evaluate (&bar, s);" << endl;
      myfile << "  " << name << "= TSIL_GetFunction (&bar,\"Suxv" <<"\");"<< endl;
    }
    if (type == "T")
    {
      string A = base.e1,B=base.e2,C=base.e3;
      myfile << "  TSIL_SetParametersST (&bar," << A << "2, " << B << "2, " << C <<"2, Q2);" << endl;
      myfile << "  TSIL_Evaluate (&bar, s);" << endl;
      myfile << "  " << name << "= -TSIL_GetFunction (&bar,\"Txuv" <<"\");"<< endl;
    }
		if (type == "TBAR")
    {
      string A = base.e1,B=base.e2,C=base.e3;
      myfile << "  TSIL_SetParametersST (&bar," << A << "2, " << B << "2, " << C <<"2, Q2);" << endl;
      myfile << "  TSIL_Evaluate (&bar, s);" << endl;
      // note minus sign accounted for in MassBuilder.m
      myfile << "  " << name << "= TSIL_GetFunction (&bar,\"TBARxuv" <<"\");"<< endl;
    }    
    if (type == "F")
    {
      string A = base.e1,B=base.e2,C=base.e3,D=base.e4,E=base.e5;
      myfile << "  TSIL_SetParameters (&bar," << A << "2, " << B << "2, " << C << "2 , " << D << "2 , " << E  << "2, Q2);" << endl;
      myfile << "  TSIL_Evaluate (&bar, s);" << endl;
      myfile << "  " << name << "= TSIL_GetFunction (&bar,\"M" <<"\");"<< endl;
    }
    if (type == "V")
    {
      string A = base.e1,B=base.e2,C=base.e3,D=base.e4;
      myfile << "  TSIL_SetParameters (&bar," << D << "2, " << C << "2, " << B << "2 , " << "1.0" << " , " << A  << "2, Q2);" << endl;
      myfile << "  TSIL_Evaluate (&bar, s);" << endl;
      myfile << "  " << name << "= -TSIL_GetFunction (&bar,\"Uzxyv" <<"\");"<< endl;
    }
    
        
    if (type == "const" || type == "A" || type == "B")
    {
			// do nothing
		}
		else
		{
			myfile << " } " << endl;
		}
    
    
  }
  
  void timestamp ( )
  {
# define TIME_SIZE 40
    
    static char time_buffer[TIME_SIZE];
    const struct std::tm *tm_ptr;
    
    std::time_t now;
    
    now = std::time ( NULL );
    tm_ptr = std::localtime ( &now );
    
    std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
    
    std::cout << time_buffer << "\n";
    
    return;
# undef TIME_SIZE
  }
  
  
  void print_diagram_info(Options options)
  {
    if ((options.particle_1!=options.particle_2))
    {
      if (options.counter_terms == true)
      {
        cout << "calculating counter-term diagram " << options.diagram << " for particle " << options.particle_1 << " to " << options.particle_2 << " in model ";
      }
      else
      {
        cout << "calculating diagram " << options.diagram << " for particle " << options.particle_1 << " to " << options.particle_2 << " in model ";
      }
    }
    else
    {
      if (options.counter_terms == true)
      {
        cout << "calculating counter-term diagram " << options.diagram << " for particle " << options.particle_1 << " in model ";
      }
      else
      {
        cout << "calculating diagram " << options.diagram << " for particle " << options.particle_1 << " in model ";
      }
    }
    
    cout << options.model << " at " << options.loop_order << "-loop order" << endl;
  }
  
  
  
}
