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








std::map <std::string, Bases > set_bases(std::vector<string> masses_input, std::vector<string> identifiers_input)
{
  vector<string> identifiers;
  vector<string> masses;
  
  // first need to make sure there are no duplicate masses, and assign a unique one character identifier to each
  // mass
    
  cout << "masses are " << endl;
  for (unsigned int i = 0; i < masses_input.size(); i++)
  {
  cout << masses_input[i] << endl;
  }
  
  
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
  
  identifiers = remove_duplicates(identifiers_input,"duplicate identifiers found");
  
  
  
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
  
  cout << "identifiers are " << endl;
  for (unsigned int i = 0; i < identifiers.size(); i++)
  {
  cout << identifiers[i] << endl;
  }
  cout << " --- " << endl;
  
  
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