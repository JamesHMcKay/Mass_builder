#ifndef DATA_H
#define DATA_H
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
struct Data
{
public:

double P, Q;
  Data (){};
Data(int argc, char* argv[]) {
double param [99];
std::string name [99]; int i=0;
if (argc==2)
{
cout << "Please enter a data file in the format ./main input.txt, using default values " << endl;
}
else
{
std::ifstream input(argv[2]);
std::string line;
while(getline(input, line)) {
if (!line.length() || line[0] == '#')
 continue;
 std::istringstream iss(line);
 iss>> name[i] >> param[i];
     i=i+1;
   }
  }
  for (int n=0;n<i+1;n++)
  {
  if (name[n]=="Q")
  {  Q =param[n];  }
  if (name[n]=="P")
  {  P =param[n];  }
}
}
};
#endif
