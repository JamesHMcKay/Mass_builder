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







}
