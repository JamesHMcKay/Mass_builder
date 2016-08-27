#include "utils.hpp"

namespace utils
{
void get_data(vector<std::string> &A,int &n,const char *filename)
{

//cout << "reading file = " << filename << endl;
n=0;
std::ifstream input(filename);
std::string line;
while(getline(input, line)) {
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
while(getline(input2, line2)) {
    if (!line2.length() || line2[0] == '#')
       continue;
    std::istringstream iss2(line2);
  
  
  iss2>> A[n];
    n=n+1;
 }

 input2.close();


}

}