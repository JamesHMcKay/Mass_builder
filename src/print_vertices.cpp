/*
 Mass Builder
 
 James McKay
 Feb 2017
 
 --- calc_counter_terms.cpp ---
 
 */

#include "calc_counter_terms.hpp"

#define RUN_ALL
//#define DEBUG
using namespace std;
using namespace utils;
string s_cwd__(getcwd(NULL,0));


void print_vertex(ofstream &myfile, std::string particle_1,std::string particle_2,std::string particle_3,Options options)
{
  std::string type;
  if (options.counter_terms){myfile<<"t12 = CreateCTTopologies[1, 1 -> 2, ExcludeTopologies -> Internal];"<<endl;
    type = "c";}
  else {myfile<<"t12 = CreateTopologies[0, 1 -> 2, ExcludeTopologies -> Internal];"<<endl;
    type = "";}
  
  myfile<<"alldiags = InsertFields[t12, {"<<particle_1<<"} -> {"<<particle_2<<","<<particle_3<<"},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions-> {" << options.restrictions << "}, Model -> \""<<s_cwd__<<"/models/"<<options.model<<"/"<<options.model<<"\"];\n";
  myfile<<"amp0 = FCFAConvert[CreateFeynAmp[alldiags], IncomingMomenta -> {Subscript[p, 1]},OutgoingMomenta -> {Subscript[p, 2], Subscript[p, 3]},UndoChiralSplittings -> True, TransversePolarizationVectors -> {p},DropSumOver -> True, List -> False, ChangeDimension -> D] //Contract\n";
  
  myfile<<"Export[\""<<s_cwd__<<"/models/"<<options.model<<"/LaTeX/vertex_"<<particle_1<< "_"<<particle_2<< "_"<<particle_3<< "_" << type <<".txt\",TeXForm[amp0]];\n"
  
  <<"Export[\""<<s_cwd__<<"/models/"<<options.model<<"/LaTeX/vertex_"<<particle_1<< "_"<<particle_2<< "_"<<particle_3<< "_" << type <<".pdf\",Paint[alldiags, Numbering -> None, ColumnsXRows -> 1]];\n"
  
  
  <<endl;
  
}


void print_vertex_tex(ofstream &myfile, std::string particle_1,std::string particle_2,std::string particle_3,Options options)
{
  std::string type;
  if (options.counter_terms){ type = "c";}
  else {type = "";}
  
  // check if the vertex has been calculated correctly
  
  
  const char* file_vertices_tmp = "models/";
  string c_file_vertices = file_vertices_tmp + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 + "_" + particle_3 + "_" + type + ".txt";
  const char *file_vertices = c_file_vertices.c_str();
  
  
  std::ifstream file(file_vertices);
  std::string str;
  std::string result;
  std::getline(file, str);
  result += str;
  // need to check if the result contains a basis integral
  bool success = true;
  
  std::string bad_list[5] = {"FCFAConvert"};
  
  for (int i = 0; i<1; i++)
  {
    if (result.find(bad_list[i]) != std::string::npos)
    {
      success = false;
    }
  }
  
  // fix \text{MajoranaSpinor} problem
  
  ifstream infile(file_vertices);

  string content = "";
  int i;

  for(i=0 ; infile.eof()!=true ; i++)
  {
    content += infile.get();
  }
  i--;
  content.erase(content.end()-1);

  infile.close();

  
  std::string from = "\\text{MajoranaSpinor}";
  std::string to = "\\psi";
  ReplaceAll(content,from, to);

  from = "^*^{\\text{Lor1}}";
  to = "^{*,\\text{Lor1}}";
  ReplaceAll(content,from, to);


  from = "\\text{SAB}";
  to = "\\sin(\\alpha+\\beta)";
  ReplaceAll(content,from, to);
  from = "\\text{CAB}";
  to = "\\cos(\\alpha+\\beta)";
  ReplaceAll(content,from, to);
  
  from = "\\text{SBA}";
  to = "\\sin(\\beta-\\alpha)";
  ReplaceAll(content,from, to);
  from = "\\text{CBA}";
  to = "\\cos(\\beta-\\alpha)";
  ReplaceAll(content,from, to);

  from = "\\text{CA}";
  to = "\\cos(\\alpha)";
  ReplaceAll(content,from, to);
  
  from = "\\text{SA}";
  to = "\\sin(\\alpha)";
  ReplaceAll(content,from, to);
  
  from = "\\text{CB}";
  to = "\\cos(\\beta)";
  ReplaceAll(content,from, to);
  
  from = "\\text{SB}";
  to = "\\sin(\\beta)";
  ReplaceAll(content,from, to);
  
  from = "\\text{TB}";
  to = "\\tan(\\beta)";
  ReplaceAll(content,from, to);
  
  from = "\\text{S2A}";
  to = "\\sin(2\\alpha)";
  ReplaceAll(content,from, to);
  
  from = "\\text{C2A}";
  to = "\\cos(2\\alpha)";
  ReplaceAll(content,from, to);
  
  from = "\\text{C2B}";
  to = "\\cos(2\\beta)";
  ReplaceAll(content,from, to);
  
  from = "\\text{S2B}";
  to = "\\sin(2\\beta)";
  ReplaceAll(content,from, to);
  
    
  from = "\\text{C2TW}";
  to = "\\cos(2\\theta_W)";
  ReplaceAll(content,from, to);
  
    
  from = "\\text{S2TW}";
  to = "\\sin(2\\theta_W)";
  ReplaceAll(content,from, to);
  
  
  from = "\\text{mw}";
  to = "m_W";
  ReplaceAll(content,from, to);
  
  from = "\\text{mz}";
  to = "m_Z";
  ReplaceAll(content,from, to);
  
  from = "\\text{MChi}";
  to = "m_{\\chi}";
  ReplaceAll(content,from, to);
  
  from = "\\text{cw}";
  to = "\\sin(\\theta_W)";
  ReplaceAll(content,from, to);
  
  from = "\\text{sw}";
  to = "\\cos(\\theta_W)";
  ReplaceAll(content,from, to);
  
  from = "\\text{STW}";
  to = "\\sin(\\theta_W)";
  ReplaceAll(content,from, to);
  
  from = "\\text{CTW}";
  to = "\\cos(\\theta_W)";
  ReplaceAll(content,from, to);
  
  from = "\\text{Lam}";
  to = "\\lambda";
  ReplaceAll(content,from, to);
  

  

  
  
  if (success)
  {
    
    myfile << "\\begin{minipage}{0.2\\textwidth}\n"
    <<"\\includegraphics[width=1\\textwidth]{"<<s_cwd__<<"/models/"<<options.model<<"/LaTeX/vertex_"<<particle_1<< "_"<<particle_2<< "_"<<particle_3<< "_" << type <<".pdf}\n"
    <<"\\end{minipage}\n"
    <<"\\hspace{2cm}\n"
    <<"\\begin{minipage}{0.7\\textwidth}\n"
    <<"\\begin{math}\n"
    << content
    //<<"\\input{"<<s_cwd__<<"/models/"<<options.model<<"/LaTeX/vertex_"<<particle_1<< "_"<<particle_2<< "_"<<particle_3<< "_" << type <<".txt}\n"
    <<"\\end{math}\n"
    <<"\\end{minipage}\n"
    <<"\n"
    <<"\\vspace{1cm}\n"
    <<"\n";
    
  }
  else
  {
    myfile << "\\begin{minipage}{0.29\\textwidth}\n"
    <<"\\small{\\lstinline{"<<particle_1<< "} $\\rightarrow$ \\lstinline{"<<particle_2<< "} $+$ \\lstinline{"<<particle_3<< "}} \n"
    <<"\\end{minipage}\n"
    <<"\\begin{minipage}{0.7\\textwidth}\n";
    if (type == "c" )
    {
      myfile<< "The requested counter-term vertex does not exist.\n";
    }
    else
    {
      myfile<< "The requested vertex does not exist.\n";
    }
    myfile<<"\\end{minipage}\n";
  }
  
  myfile<<"\n"
  <<"\\vspace{1cm}\n"
  <<"\n";
  
}



void print_vertices(Options options)
{
  
  
  
  
  string s_cwd__(getcwd(NULL,0));
  
  string model = options.model;
  
  
  // read in a list of 3-particle vertices from file
  
  vector<string> particle_1,particle_2,particle_3;
  int np;
  
  const char* file_vertices_tmp = "models/";
  string c_file_vertices = file_vertices_tmp + model + "/vertices.txt";
  const char *file_vertices = c_file_vertices.c_str();
  
  
  get_data(particle_1,particle_2,particle_3,np, file_vertices);
  
  
  // create Mathematica script to print a pdf for each vertex and record the amplitude in a text file
  
  ofstream myfile;
  myfile.open("output/print_vertices.m");
  
  utils::print_math_header(myfile);
  utils::assign_FCGV(myfile,options);
  utils::assign_variables(myfile,options);
  
  // loop over subroutine for all particles in list
  
  for (int i = 0; i < np; i++)
  {
    print_vertex(myfile,particle_1[i],particle_2[i],particle_3[i],options);
  }
  
  // run Mathematica script
  
  myfile.close();
  
#ifdef RUN_ALL
  system("chmod +x output/print_vertices.m ");
  if (options.verbose) system("./output/print_vertices.m");
  else system("./output/print_vertices.m  >/dev/null");
#endif
  
  
  
  
  // generate LaTeX file
  
  const char* file_tex_tmp = "models/";
  string c_file_tex = file_tex_tmp + model + "/LaTeX/vertices.tex";
  const char *file_tex = c_file_tex.c_str();
  
  ofstream texfile;
  texfile.open(file_tex);
  
  texfile << "\\documentclass[11pt]{amsart}\n"
  << "\\usepackage{geometry}\n"
  << "\\geometry{letterpaper}\n"
  
  << "\\usepackage[parfill]{parskip}\n"
  << "\\usepackage{graphicx}\n"
  << "\\usepackage{amssymb}\n"
  << "\\usepackage{epstopdf}\n"
  <<"\\usepackage{listings}\n"
  <<"\\usepackage{breqn}\n"
  << "\\DeclareGraphicsRule{.tif}{png}{.png}{`convert #1 `dirname #1`/`basename #1 .tif`.png}\n"
  
  << "\\title{\\lstinline{Vertices in the "<<options.model<<"}}\n"
  << "\\author{File generated by \\textsf{Mass Builder} on \\today}\n"
  << "%\\date{\\today}\n"
  
  << "\\begin{document}\n"
  
  << "\\maketitle\n";
  
  
  
  for (int i = 0; i < np; i++)
  {
    print_vertex_tex(texfile,particle_1[i],particle_2[i],particle_3[i],options);
  }
  
  
  texfile<<"\\end{document}\n";
  
  texfile.close();
  
}