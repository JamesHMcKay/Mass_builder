/* 

The functions defined in this file add the ability to significantly reduce the number TSIL_Evaluate calls required during
a self energy calculation by taking advantage of the symmetries between basis integrals in an entirely
generic way

Sep 2016

*/

#include "write_tsil_ini.hpp"


void Print_dotsil::print_to_file(ofstream &myfile)
{
  sort_integrals();
  
}


void Print_dotsil::sort_integrals()
{
// start with

// find all the TSIL_EVALUATE calls that will satisfy a particular basis integral




}



void Print_dotsil::get_poss_eval(Bases base)
{

// for each type of integral find all the evaluate_tsil statements that would provide it

if (base.type == "V")
{
// for V(u,x,z,y) the corresponding TSIL type is -U(x,y,z,u)




}





}