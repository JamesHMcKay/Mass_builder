#ifndef DATA_H
#define DATA_H
#include "options.hpp"
using namespace std;
struct Data
{
public:
double sw;
double cw;
double sw2;
double cw2;
double S2TW;
double g2;
double g1;
double CTW;
double C2TW;
double MChi;
double ma;
double mw;
double mz;
double mg;
double mh;
double MassVZ;
double MassVWp;
double MassAh;
double MassHp;
double SE_F5;
double SE_F6;
double P, Q;
  Data (){};
Data(Options options) {
double param [99];
std::string name [99]; int i=0;
std::ifstream input(options.input_list);
std::string line;
while(getline(input, line)) {
if (!line.length() || line[0] == '#')
 continue;
 std::istringstream iss(line);
 iss>> name[i] >> param[i];
     i=i+1;
   }
  for (int n=0;n<i+1;n++)
  {
  if (name[n]=="sw")
  {
  sw=param[n];
  }
  if (name[n]=="cw")
  {
  cw=param[n];
  }
  if (name[n]=="sw2")
  {
  sw2=param[n];
  }
  if (name[n]=="cw2")
  {
  cw2=param[n];
  }
  if (name[n]=="S2TW")
  {
  S2TW=param[n];
  }
  if (name[n]=="g2")
  {
  g2=param[n];
  }
  if (name[n]=="g1")
  {
  g1=param[n];
  }
  if (name[n]=="CTW")
  {
  CTW=param[n];
  }
  if (name[n]=="C2TW")
  {
  C2TW=param[n];
  }
  if (name[n]=="MChi")
  {  MChi=param[n];  }
  if (name[n]=="ma")
  {  ma=param[n];  }
  if (name[n]=="mw")
  {  mw=param[n];  }
  if (name[n]=="mz")
  {  mz=param[n];  }
  if (name[n]=="mg")
  {  mg=param[n];  }
  if (name[n]=="mh")
  {  mh=param[n];  }
  if (name[n]=="MassVZ")
  {  MassVZ=param[n];  }
  if (name[n]=="MassVWp")
  {  MassVWp=param[n];  }
  if (name[n]=="MassAh")
  {  MassAh=param[n];  }
  if (name[n]=="MassHp")
  {  MassHp=param[n];  }
  if (name[n]=="Q")
  {  Q =param[n];  }
  if (name[n]=="P")
  {  P =param[n];  }
}
}
};
#endif
