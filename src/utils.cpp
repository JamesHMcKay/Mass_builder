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
  
  int get_loop_order(string type)
  {
    int loop_order;
    string loop_order_str = utils::char_to_string(type[0]);
    stringstream convert(loop_order_str);
    convert >> loop_order;
    return loop_order;
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
  
  
  
  
  void user_input_guide()
  {
    cout << " Welcome to Mass Builder \n"
    << " ---------------------------------------------------------------------------------------------------- \n"
    << " to use call ./mass_builder <flags> where required flags for each mode are given below:\n"
    << "                                                                                                      \n"
    << " -a -m <model>                                           -  compute all diagrams in models/<model>/diagrams.txt\n"
    << " -a -m <model>                               -i <file>   -  compute all diagrams in listed in file\n"
    << " -a -m <model>  -p <particle>  -d <diagram>              -  compute specific diagram\n"
    << " -g -m <model>                                           -  generate code for available diagrams\n"
    << " -g -m <model>                               -i <file>   -  generate code for diagrams listed in file\n"
    << " -f -m <model>  -p <particle>                            -  draw all FeynArts diagrams for particle\n"
    << " -e                                          -i <file>   -  evaluate self energy with values for masses and couplings in file\n"
    << " ---------------------------------------------------------------------------------------------------- " <<endl;
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
  
  
  void print_math_header(ofstream &file)
  {
    /*MATH_PATH */  file<< "#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script"<<endl;
    file << "(* ---------------------------------------------------- *)\n"
    << "(* This file has been automatically generated by Mass Builder, on "<< now->tm_mday << '-'
    << (now->tm_mon + 1) << '-'<< (now->tm_year + 1900) <<"*)\n"
    << "(* ---------------------------------------------------- *)\n"
    << "$LoadTARCER = True;\n"
    <<"$LoadFeynArts = True;\n"
    << "<< FeynCalc/FeynCalc.m\n"
    <<"dm[mu_] := DiracMatrix[mu, Dimension -> D]\n"
    <<"dm[5] := DiracMatrix[5]\n"
    <<"ds[p_] := DiracSlash[p]\n"
    <<"SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True];\n"
    <<"SetOptions[DiracTrace, DiracTraceEvaluate -> True];\n"
    <<"$GenericMixing = True;\n";
  }
 
void print_math_body(ofstream &file,Options options,string cwd,std::vector<std::string> masses)
  {
    int loop_order = options.loop_order;
    string particle_1 = options.particle_1;
    string particle_2 = options.particle_2;
    string diagram = options.diagram;
    string model = options.model;
    if (options.counter_terms)
    {
      file<<"t12 = CreateCTTopologies["<< loop_order<<", 1 -> 1, ExcludeTopologies -> Internal];"<<endl;
    }
    else
    {
      file<<"t12 = CreateTopologies["<< loop_order<<", 1 -> 1, ExcludeTopologies -> Internal];"<<endl;
    }
    
    file <<"alldiags = InsertFields[t12, {"<<particle_1<<"} -> {"<<particle_2<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Model -> \""<<cwd<<"/models/"<<model<<"/"<<model<<"\"];\n"
    <<"subdiags0 =   DiagramExtract[alldiags, "<<diagram<<"]\n"
    <<"amp0 = FCFAConvert[CreateFeynAmp[subdiags0], IncomingMomenta -> {p}, OutgoingMomenta -> {p}, LoopMomenta -> {k1, k2} ,UndoChiralSplittings -> True,TransversePolarizationVectors -> {p},DropSumOver -> True, List -> False,ChangeDimension -> D] // Contract\n";
    for (unsigned int i=0; i < masses.size();i++)
    {
      file<<"amp0 = amp0 /. MajoranaSpinor[Momentum[p, D], "<<masses[i]<<"] -> 1 /.Spinor[Momentum[p, D], "<<masses[i]<<", 1] -> 1;"<<endl;
    }
    file<<"SetOptions[Eps, Dimension -> D];\n";
    
    if ( (loop_order == 2) && (!options.counter_terms) )
    {
      file<<"fullamp0 = (amp0) // DiracSimplify // FCMultiLoopTID[#, {k1, k2}] & //DiracSimplify;\n"
      <<"tfiamp0 = fullamp0 // ToTFI[#, k1, k2, p] & // ChangeDimension[#, D] &;\n";
    }
    
    else if ( (loop_order == 1) && (options.counter_terms) )
    {
      file<<" fullamp0 = (amp0) // DiracSimplify;\n"
      <<"tfiamp0 = fullamp0 // ChangeDimension[#, D] &;\n";
    }
    else
    {
      file<<" fullamp0 = (amp0) // DiracSimplify // TID[#, k1] & // DiracSimplify;\n"
      <<"tfiamp0 = fullamp0 // ToTFI[#, k1, p] & // ChangeDimension[#, D] &;\n";
    }
    
    for (unsigned int i=0; i < masses.size();i++)
    {
      file<<"tfiamp0 = tfiamp0 /. MajoranaSpinor[Momentum[p,D], "<<masses[i]<<"] -> 1 /.Spinor[Momentum[p,D], "<<masses[i]<<", 1] -> 1;"<<endl;
      file<<"tfiamp0 = tfiamp0 /. MajoranaSpinor[p, "<<masses[i]<<"] -> 1 /.Spinor[Momentum[p], "<<masses[i]<<", 1] -> 1;"<<endl;
    }
    
  }
  
  
  bool check_done_quiet()
  {
    std::ifstream file("output/result.txt");
    std::string str;
    std::string result;
    std::getline(file, str);
    result += str;
    // need to check if the result contains a basis integral
    bool success = true;
    
    std::string bad_list[5] = {"TFI","TBI","TVI","TAI","TJI"};
    
    for (int i = 0; i<5; i++)
    {
      if (result.find(bad_list[i]) != std::string::npos)
      {
        success = false;
      }
    }


    return success;
  }
  
  bool check_done(string &result)
  {
    std::ifstream file("output/result.txt");
    std::string str;
    std::getline(file, str);
    result += str;
    // need to check if the result contains a basis integral
    bool success = true;
    
    std::string bad_list[5] = {"TFI","TBI","TVI","TAI","TJI"};
    
    for (int i = 0; i<5; i++)
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
  
  
  std::string part_name_simple(std::string particle_name_full)
  {
    stringstream _part_1,_part_2;
    string part_1,part_2;
    
    _part_1 << particle_name_full[0];
    _part_1 >> part_1;
    _part_2 << particle_name_full[2];
    _part_2 >> part_2;
    
    return part_1+part_2;
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
  
  
  void print_base(ofstream &myfile, Bases base, string id, string SEn, string D)
  {
    string type = base.type;
    string momentum = "Pair[Momentum[p],Momentum[p]]";
    if (D=="D")
    {
      momentum = "Pair[Momentum[p,D],Momentum[p,D]]";
    }
    
    if (SEn != "diff")
    {
      if (type == "F")
      {
        myfile << id << " = " << "TFI["<<D<<", " << momentum << ", {{1, " << base.e1  << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}, {1, " << base.e4 << "}, {1, " << base.e5 << "}}];" << endl;
      }
      if (type == "A")
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
      if (type == "J")
      {
        myfile << id << " = " << "TJI["<<D<<", " << momentum << ", {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      }
      if (type == "K")
      {
        myfile << id << " = " << "TJI["<<D<<", 0, {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      }
    }
    
    myfile << "C"<< id << " = Coefficient["<<SEn<<", " << id << ", 1];" << endl;
    myfile << "C"<< id << "2 = Coefficient["<<SEn<<", " << id << ", 2];" << endl; // check if the squared integral exists
  }
  
  
  
  
  
  
  // print the basis integrals out in Mathematica notation
  void print_math_basis(std::map<std::string, Bases> base_map, ofstream &myfile, string target, string D)
  {
    vector<string> bases_names = extract_keys(base_map);
    for (unsigned int i = 0; i < bases_names.size();i++)
    {
      Bases base_temp;
      base_temp = base_map[bases_names[i]];
      print_base(myfile, base_temp, bases_names[i], target,D);
    }
  }
  
  
  
  
  
  
  void print_base_product(ofstream &myfile,Bases base_1,Bases base_2,string SEn, string D)
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
    if (type1=="K") myfile << "TJI["<<D<<", 0, {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}}];";
    if (type1=="V") myfile << "TVI["<<D<<", " << momentum << ", {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "}}]";
    //if (type1=="F") myfile << "TFI["<<D<<", Pair[Momentum[p,D],Momentum[p,D]], {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "},{1, " << base_1.e5 << "}}]";
    myfile << " * ";
    if (type2=="A") myfile << "TAI["<<D<<", 0, {{1, " << base_2.e1 << "}}];";
    if (type2=="B") myfile << "TBI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}}];";
    if (type2=="J") myfile << "TJI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];";
    if (type2=="T") myfile << "TJI["<<D<<", " << momentum << ", {{2, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];";
    if (type2=="K") myfile << "TJI["<<D<<", 0, {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}}];;";
    if (type2=="V") myfile << "TVI["<<D<<", " << momentum << ", {{1, " << base_2.e1 << "}, {1, " << base_2.e2 << "}, {1, " << base_2.e3 << "}, {1, " << base_2.e4 << "}}];";
    //if (type2=="F") myfile << "TFI["<<D<<", Pair[Momentum[p,D],Momentum[p,D]], {{1, " << base_1.e1 << "}, {1, " << base_1.e2 << "}, {1, " << base_1.e3 << "}, {1, " << base_1.e4 << "},{1, " << base_1.e5 << "}}]";
    myfile << endl;
    if (base_1.short_name==base_2.short_name) myfile << "C"<< base_1.short_name << base_2.short_name << " = Coefficient["<<SEn<<","<< base_1.short_name << ", 2];" << endl;
    else myfile << "C"<< base_1.short_name << base_2.short_name << " = - (1/2)* Coefficient["<<SEn<<","<< base_1.short_name << base_2.short_name << ", 1];" << endl;
  end:;
  }
  
  
  void print_math_products(std::map<std::string, Bases> base_map, ofstream &myfile, string target, string D)
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
        
        print_base_product(myfile,base_1,base_2,target,D);
        
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
    if (type == "const")
    {
      myfile << "  const" << " = 1.0L"<<endl;
    }
    if (type == "A")
    {
      myfile << "  " <<  name << " = -i*TSIL_A_ ("<<base.e1<<"2 , Q2);"<<endl;
    }
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
    
  }
  
  
  
  // expressing basis integrals as finite + divergent parts
  // x -> x^2 unless in a TXI function call
  // TSIL -> TARCER
  // A(x)       =   ia TAI(x,y)
  // B(x,y)     = - ia TBI(x,y)
  // I(x,y,z)   =   a^2 TJI[0,{1,1,1}] = a^2 K(x,y,z)
  // S(x,y,z)   =   a^2 TJI[s,{1,1,1}] = a^2 J(x,y,z)
  // T(x,y,z)   = - a^2 TJI[s,{2,1,1}] = a^2 T(x,y,z)
  // U(x,y,z,u) = - a^2 TVI[s,{1,1,1,1},{u,x,z,y}] = a^2 V(u,x,z,y)
  
  
  void finite_A(ofstream &myfile, string x, string id)
  {
    myfile << id << " = " <<  id << "4 +  ( ";
    
    myfile << " I*" << x <<"^2 /epsilon";
    myfile << ")";
    
    if (add_epsilon_terms)
    {
      myfile << " - I * epsilon * Ae[ " << x << " ] " << endl;
    }
    else
    {
      myfile << endl;
    }
    
  }
  
  void finite_B(ofstream &myfile, string x, string y, string id)
  {
    myfile << id << " = " <<  id << "4 +  ( ";
    
    myfile << " I/epsilon";
    myfile << ")";
    
    if (add_epsilon_terms)
    {
      myfile << " + I * epsilon * Be[" << x << " , " << y << "]" << endl;
    }
    else
    {
      myfile << endl;
    }
    
  }
  
  void finite_J(ofstream &myfile, string x, string y, string z, string id)
  {
    myfile << id << " = " <<  id << "4 + ( ";
    
    myfile << " - (" << x << "^2 + " << y << "^2 + " << z << "^2 )/ (2* epsilon^2)";
    myfile << " + (- ("<<x<<"^2 +"<<y<<"^2 +"<<z<<"^2)/2 + Pair[Momentum[p],Momentum[p]]/4 )/epsilon";
    myfile << " + ( I*TAI[4, 0, {1, " << x << "}] + I*TAI[4, 0, {1, " << y << "}] + I*TAI[4, 0, {1, " << z << "}] )/epsilon";
    myfile << " ) ";
    
    if (add_epsilon_terms)
    {
      myfile << " + Ae[ " << x << " ] + " << " Ae[ " << y << " ] + Ae[ " << z << " ] " << endl;
    }
    else
    {
      myfile << endl;
    }
    
  }
  
  void finite_T(ofstream &myfile, string x, string y, string z, string id)
  {
    myfile << id << " = " <<  id << "4 + ( ";
    
    myfile << " -1/(2*epsilon^2) + 1/(2*epsilon)";
    myfile << " + ((I*TAI[4, 0, {1, "<<x<<"}])/"<<x<<"^2 )/epsilon"; // was the x^2 not being ^ a typo? very late change here
    myfile << " - ((I*TAI[4, 0, {1, "<<x<<"}])/"<<x<<"^2 )";  // was the x^2 not being ^ a typo? very late change here
    myfile << " ) ";
    
    if (add_epsilon_terms)
    {
      myfile << " + Ae[ " << x << " ]/"<< x <<"^2" << endl;
    }
    else
    {
      myfile << endl;
    }
    
  }
  
  void finite_K(ofstream &myfile, string x, string y, string z, string id)
  {
    myfile << id << " = " <<  id << "4 +  ( ";
    
    myfile << " - (" << x << "^2 + " << y << "^2 + " << z << "^2 )/ (2* epsilon^2)";
    myfile << " - (" << x << "^2 + " << y << "^2 + " << z << "^2 )/  ( 2*epsilon)";
    myfile << "+  (  I*TAI[4, 0, {1, "<<x<<"}]";
    myfile << "  + I*TAI[4, 0, {1, "<<y<<"}]";
    myfile << "  + I*TAI[4, 0, {1, "<<z<<"}]";
    myfile << "    )/epsilon ";
    myfile << " ) ";
    
    if (add_epsilon_terms)
    {
      myfile << " + Ae[ " << x << " ] + " << " Ae[ " << y << " ] + Ae[ " << z << " ] " << endl;
    }
    else
    {
      myfile << endl;
    }
    
    
  }
  
  
  void finite_V(ofstream &myfile, string x, string y, string z, string u, string id)
  {
    myfile << id << " = " <<  id << "4 + ( ";
    
    myfile << "-1/(2*epsilon^2) - 1/(2*epsilon)";
    myfile << " - ( - I * TBI[4, Pair[Momentum[p],Momentum[p]], {{1, " << y << "}, {1, " << u << "}}] ) /epsilon";
    myfile << " ) ";
    
    if (add_epsilon_terms)
    {
      myfile << " - Be[" << y << " , " << u << "]" << endl;
    }
    else
    {
      myfile << endl;
    }
    
    
  }

  
  void print_finite_base(ofstream &myfile, Bases base, string id)
  {
    string type = base.type;
    
    if (type == "F")
    {
      myfile << id << "4 = " << "TFI[4, Pair[Momentum[p],Momentum[p]], {{1, " << base.e1  << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}, {1, " << base.e4 << "}, {1, " << base.e5 << "}}];" << endl;
      myfile << id << " = " << id << "4 " << endl;
    }
    if (type == "A")
    {
      myfile << id << "4 = " <<  "TAI[4, 0, {1, " << base.e1 << "}];" << endl;
      finite_A(myfile,base.e1,id);
    }
    if (type == "B")
    {
      myfile << id << "4 = " << "TBI[4, Pair[Momentum[p],Momentum[p]], {{1, " << base.e1 << "}, {1, " << base.e2 << "}}];" << endl;
      finite_B(myfile,base.e1,base.e2,id);
    }
    if (type == "V")
    {
      myfile << id << "4 = " << "TVI[4, Pair[Momentum[p],Momentum[p]], {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}, {1, " << base.e4 << "}}];" << endl;
      finite_V(myfile, base.e1,base.e2,base.e3,base.e4,id);
    }
    if (type == "T")
    {
      myfile << id << "4 = " << "TJI[4, Pair[Momentum[p],Momentum[p]], {{2, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      finite_T(myfile,base.e1,base.e2,base.e3,id);
    }
    if (type == "J")
    {
      myfile << id << "4 = " << "TJI[4, Pair[Momentum[p],Momentum[p]], {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      finite_J(myfile,base.e1,base.e2,base.e3,id);
    }
    if (type == "K")
    {
      myfile << id << "4 = " << "TJI[4, 0, {{1, " << base.e1 << "}, {1, " << base.e2 << "}, {1, " << base.e3 << "}}];" << endl;
      finite_K(myfile,base.e1,base.e2,base.e3,id);
    }
  }
  
  
  // print the basis integrals out in Mathematica notation
  void print_finite_basis(std::map<std::string, Bases> base_map, ofstream &myfile)
  {
    vector<string> bases_names = extract_keys(base_map);
    for (unsigned int i = 0; i < bases_names.size();i++)
    {
      Bases base_temp;
      base_temp = base_map[bases_names[i]];
      print_finite_base(myfile, base_temp, bases_names[i]);
    }
  }
  
  
  
  
}
