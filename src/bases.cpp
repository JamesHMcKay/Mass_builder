#include "bases.hpp"
#include "utils.hpp"
#include <exception>
#include <sstream>
#include <cstdlib>
#include <string>

// compute all possible bases integral permutations given a list of masses return as a map
// consisting of short name and bases object
using namespace std;
using namespace utils;

// simple helper function to search for duplicates and return vector of unique elements


void set_id(std::vector<string> &masses_input, std::vector<string> &identifiers_input)
{
  vector<string> identifiers;
  vector<string> masses;
  
  // first need to make sure there are no duplicate masses, and assign a unique one character identifier to each
  // mass
  #ifdef DEBUG
  cout << "masses are " << endl;
  for (unsigned int i = 0; i < masses_input.size(); i++)
  {
  cout << masses_input[i] << endl;
  }
  #endif
  
  if (identifiers_input.size()==0)  // need to create unique identifiers
  {
  
  masses = remove_duplicates(masses_input,"Duplicate input masses duplicates are removed but please check this was not a case of a mislabelled mass and that all masses are accounted for.");
  
  cout << " no identifiers provided so trying to make them up from given masses " << endl;
  // now duplicates have been removed we will try and assign a one character name to each mass
  
  // attempt to make list now, rule: if one char long use that, if two use the second
  
  vector<int> lengths = utils::find_string_lengths(masses);
  
  // create new list of identifiers
  for (unsigned int i = 0; i < masses.size(); i++)
  {
  string mass = masses[i];
  if ((lengths[i]) == 1){ identifiers.push_back( char_to_string(mass[0]) ); }
  if ((lengths[i]) == 2){ identifiers.push_back( char_to_string(mass[1]) ); }
  if ((lengths[i]) > 2){ identifiers.push_back( char_to_string(mass[1]) + char_to_string(mass[2]) );
  cout << "WARNING: mass name " << mass <<" is long and will result in ugly code with long basis integral identifiers" << endl; }
  }
  
  // check for duplicates in new list
  if ( !(remove_duplicates(identifiers).size()) == (identifiers.size()) )
  {
  // if the above hasn't worked just throw an error, can do something fancy later
  cout << "FATAL ERROR, mass names are too long and too similar please provide unique identifiers or relabel, if two masses are intended to be the same then give them the same identifier, this is the only way for Mass_builder to know it's safe to assign these the same mass";
  exit(1);
  }
  
  }
  else  // here we deal with case of identifiers being provided, check if list is correct size and no duplicates
  {
  
  identifiers = identifiers_input;
  
  if (remove_duplicates(identifiers_input).size()!=masses_input.size())
  {
  cout << "ERROR: duplicate names in the identifiers list, Mass Builder is not yet set up to handle this (work in progress)" << endl;
  exit(1);
  }
  
  
  vector<int> lengths = utils::find_string_lengths(identifiers);
  
  // create new list of identifiers
  for (unsigned int i = 0; i < identifiers.size(); i++)
  {
  if (lengths[i] > 1)
  {
  cout << "WARNING: identifier name " << identifiers[i] << " is long and will result in ugly code with long basis integral identifiers" << endl; }
  }
  
  masses = masses_input;
  
  
  }

  
  vector<string> id = identifiers;
  #ifdef DEBUG
  cout << "identifiers are " << endl;
  for (unsigned int i = 0; i < identifiers.size(); i++)
  {
  cout << identifiers[i] << endl;
  }
  cout << " --- " << endl;
  #endif
  
  identifiers_input = identifiers;
  masses_input = masses;

}








std::map <std::string, Bases > set_bases(std::vector<string> masses, std::vector<string> &identifiers_input)
{
  set_id( masses, identifiers_input);
  vector<string> id = identifiers_input;
  
  
  std::map <std::string, Bases > bases_map;
  int n = id.size();
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  Bases a("A",masses[i1]);
  bases_map["A"+id[i1]] = a ;
  for (int i2 = 0; i2 < n ; i2++)
  {
  for (int i3 = 0; i3 < n ; i3++)
  {
  for (int i4 = 0; i4 < n ; i4++)
  {
  Bases v("V",masses[i1],masses[i2],masses[i3],masses[i4]);
  bases_map["V"+id[i1]+id[i2]+id[i3]+id[i4]]= v;
  for (int i5 = 0; i5 < n ; i5++)
  {
  Bases f("F",masses[i1],masses[i2],masses[i3],masses[i4],masses[i5]);
  bases_map["F"+id[i1]+id[i2]+id[i3]+id[i4]+id[i5]]= f;
  }}}}}
  
  
  
   //deal with symmetries appropriately
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  Bases b("B",masses[i1],masses[i2]);
  bases_map["B"+id[i1]+id[i2]]= b;
  
  }}
  
  
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = i1; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  Bases j("J",masses[i1],masses[i2],masses[i3]);
  bases_map["J"+id[i1]+id[i2]+id[i3]]= j;
  Bases k("K",masses[i1],masses[i2],masses[i3]);
  bases_map["K"+id[i1]+id[i2]+id[i3]]= k;
  
  }}}
  
    
  for (int i1 = 0; i1 < n ; i1++)
  {
  for (int i2 = 0; i2 < n ; i2++)
  {
  for (int i3 = i2; i3 < n ; i3++)
  {
  Bases t("T",masses[i1],masses[i2],masses[i3]);
  bases_map["T"+id[i1]+id[i2]+id[i3]]= t;
  }
  }
  }
  
  

  return bases_map;



}



std::map <std::string, Bases > remove_zeros(std::map <std::string, Bases > base_map, std::vector<std::string> bases_names)
{

std::map <std::string, Bases > new_base_map;

for (unsigned int i = 0; i < bases_names.size();i++)
{

string coeff = base_map[bases_names[i]].coefficient;

if ( char_to_string(coeff[2]) != "0") // the character choosen depends on the amount of white space used
{
new_base_map[bases_names[i]] = base_map[bases_names[i]];

}
}

return new_base_map;
}


// reformat the coefficient to change Mathematica expressions into readable C++ input
void format_coeff(std::map <std::string, Bases > &base_map, std::vector<std::string> bases_names,std::vector<std::string> masses, std::vector<std::string> id)
{
int nb = bases_names.size();
int nm = masses.size();

// deal with TAI and TBI objects (should extend to deal with all possible)
string from="",to="";

for (int k = 0; k<nb ; k++)
{
string coefficient = base_map[bases_names[k]].coefficient;




for (int i = 0; i<nm; i++)
{
for (int j = 0; j<nm; j++)
{
from = "TBI(4,Power(p,2),List(List(1,"+masses[i]+"),List(1," + masses[j] + ")))";  // TODO make short name for all bases object and for the B object include both permutations, then this solves this problem here
to = "B"+id[j]+id[i]; // choice of i and j here is important TODO

ReplaceAll(coefficient,from, to);
}

from = "TAI(4,0,List(List(1,"+masses[i]+")))";
to = "A"+id[i];
ReplaceAll(coefficient,from, to);

from = "MajoranaSpinor(p,"+masses[i]+")";
to = "1.0";
ReplaceAll(coefficient,from, to);

from = "Spinor(Momentum(p),"+masses[i]+",1)";
to = "1.0";
ReplaceAll(coefficient,from, to);





}

base_map[bases_names[k]].coefficient = coefficient;
}

}




std::map <std::string, Bases > products_container(vector<string> bases_names)
{
std::map <std::string, Bases > prod_map;
int n = bases_names.size();

for (int i = 0; i < n ; i++)
{
for (int j = 0; j < n ; j++)
{
string name = bases_names[i]+bases_names[j];
Bases base;
base.short_name = name;
base.e1 = bases_names[i];
base.e2 = bases_names[j];

prod_map[name] = base;
}
}

return prod_map;
}







