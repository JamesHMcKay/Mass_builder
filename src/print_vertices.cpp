/*
 Mass Builder
 
 James McKay
 May 2017
 
 --- print_vertices.cpp ---
 
 */

#include "utils.hpp"
#include "templates.hpp"
#include "compute_amp.hpp"

#define RUN_ALL
//#define DEBUG
using namespace std;
using namespace utils;


void Compute_amp::print_3_vertex(std::string &input, std::string particle_1,std::string particle_2,std::string particle_3,Options options)
{
  std::string type;
  if (options.counter_terms){input += "t12 = CreateCTTopologies[1, 1 -> 2, ExcludeTopologies -> Internal];";
    type = "c";}
  else {input += "t12 = CreateTopologies[0, 1 -> 2, ExcludeTopologies -> Internal];";
    type = "";}
  
  input += "alldiags = InsertFields[t12, {" + particle_1 + "} -> {" + particle_2 + "," + particle_3 + "},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions-> {"  +  options.restrictions  +  "}, Model -> \"" + get_cwd() + "/models/" + options.model + "/" + options.model + "\"];";
  input += "amp0 = FCFAConvert[CreateFeynAmp[alldiags], IncomingMomenta -> {Subscript[p, 1]},OutgoingMomenta -> {Subscript[p, 2], Subscript[p, 3]},UndoChiralSplittings -> True, TransversePolarizationVectors -> {p},DropSumOver -> True, List -> False, ChangeDimension -> D] //Contract;";
  
  input += "Export[\"" + get_cwd() + "/models/" + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 +  "_" + particle_3 +  "_"  +  type  + ".txt\",TeXForm[amp0]];";
  
  input += "Export[\"" + get_cwd() + "/models/" + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 +  "_" + particle_3 +  "_"  +  type  + ".pdf\",Paint[alldiags, Numbering -> None, ColumnsXRows -> 1]];";
  
}

void Compute_amp::print_4_vertex(std::string &input, std::string particle_1,std::string particle_2,std::string particle_3,std::string particle_4,Options options)
{
  std::string type;
  if (options.counter_terms){input += "t12 = CreateCTTopologies[1, 2 -> 2, ExcludeTopologies -> Internal];";
    type = "c";}
  else {input += "t12 = CreateTopologies[0, 2 -> 2, ExcludeTopologies -> Internal];";
    type = "";}
  
  input += "alldiags = InsertFields[t12, {" + particle_1 + "," + particle_2 + "} -> {" + particle_3 + "," + particle_4 + "},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions-> {"  +  options.restrictions  +  "}, Model -> \"" + get_cwd() + "/models/" + options.model + "/" + options.model + "\"];";
  input += "amp0 = FCFAConvert[CreateFeynAmp[alldiags], IncomingMomenta -> {Subscript[p, 1],Subscript[p,2]},OutgoingMomenta -> {Subscript[p, 3], Subscript[p, 4]},UndoChiralSplittings -> True, TransversePolarizationVectors -> {p},DropSumOver -> True, List -> False, ChangeDimension -> D] //Contract;";
  
  input += "Export[\"" + get_cwd() + "/models/" + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 +  "_" + particle_3 +  "_" + particle_4  +  type  + ".txt\",TeXForm[amp0]];";
  
  input += "Export[\"" + get_cwd() + "/models/" + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 +  "_" + particle_3 +  "_" + particle_4  +  type  + ".pdf\",Paint[alldiags, Numbering -> None, ColumnsXRows -> 1]];";
  
}


void Compute_amp::print_vertex_tex(ofstream &myfile, std::string particle_1,std::string particle_2,std::string particle_3,std::string particle_4,Options options,int number)
{
  std::string type;
  if (options.counter_terms){ type = "c";}
  else {type = "";}
  
  // check if the vertex has been calculated correctly
  
  
  const char* file_vertices_tmp = "models/";
  string c_file_vertices = file_vertices_tmp + options.model + "/LaTeX/vertex_" + particle_1 +  "_" + particle_2 + "_" + particle_3 + "_" + particle_4 + type + ".txt";
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
  
    
  from = "\\text{g1}";
  to = "g_1";
  ReplaceAll(content,from, to);
  
  from = "\\text{g2}";
  to = "g_2";
  ReplaceAll(content,from, to);

  

  
  
  if (success)
  {
    myfile  <<  number + 1 << ".\n"
    << "\\begin{minipage}{0.2\\textwidth}\n"
    <<"\\includegraphics[width=1\\textwidth]{"<< get_cwd() <<"/models/"<<options.model<<"/LaTeX/vertex_"<<particle_1<< "_"<<particle_2<< "_"<<particle_3<< "_"<<particle_4 << type <<".pdf}\n"
    <<"\\end{minipage}\n"
    <<"\\hspace{2cm}\n"
    <<"\\begin{minipage}{0.7\\textwidth}\n"
    <<"\\begin{math}\n"
    << content
    <<"\\end{math}\n"
    <<"\\end{minipage}\n"
    <<"\n"
    <<"\\vspace{1cm}\n"
    <<"\n";
    
  }
  else
  {
    myfile << number + 1 << ".\n"
    <<"\\begin{minipage}{0.29\\textwidth}\n";
    if (particle_4=="")
    {
      myfile<<"\\small{\\lstinline{"<<particle_1<< "} $\\rightarrow$ \\lstinline{"<<particle_2<< "} $+$ \\lstinline{"<<particle_3<< "}} \n";
    }
    else
    {
      myfile<<"\\small{\\lstinline{"<<particle_1<< "} $+$ \\lstinline{"<<particle_2<< "} $\\rightarrow$ \\lstinline{"<<particle_3<< "} $+$ \\lstinline{"<<particle_4<< "}} \n";
    }
    myfile<<"\\end{minipage}\n"
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



void Compute_amp::print_vertices()
{
  string model = options.model;
  
  
  // read in a list of 3-particle vertices from file
  
  vector<string> particle_1,particle_2,particle_3,particle_4;
  int np;
  
  const char* file_vertices_tmp = "models/";
  string c_file_vertices = file_vertices_tmp + model + "/vertices.txt";
  const char *file_vertices = c_file_vertices.c_str();
  
  
  get_data(particle_1, particle_2, particle_3, particle_4, np, file_vertices);
  
  
  // create Mathematica script to print a pdf for each vertex and record the amplitude in a text file
  
  create_wstp_link();
  load_libraries();
  WSNewPacket(link);
  std::string input;
  
  templates::print_math_header(input);
  utils::assign_FCGV(input,options);
  utils::assign_variables(input,options);
  
  send_to_math(input);
  
  // loop over subroutine for all particles in list
  
  for (int i = 0; i < np; i++)
  {
    if (particle_4[i] == "")
    {
      print_3_vertex(input,particle_1[i],particle_2[i],particle_3[i],options);
    }
    else
    {
      print_4_vertex(input,particle_1[i],particle_2[i],particle_3[i],particle_4[i],options);
    }
  }
  
  // run Mathematica script
  
  send_to_math(input);
  
  input += "Quit[]";
  send_to_math(input);
  
  WSClose(link);
  cout << "WSTP link closed successfully" << endl;
  
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
    print_vertex_tex(texfile,particle_1[i],particle_2[i],particle_3[i],particle_4[i],options,i);
  }
  
  
  texfile<<"\\end{document}\n";
  
  texfile.close();
  
}
