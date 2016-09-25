/* 

The functions defined in this file add the ability to significantly reduce the number TSIL_Evaluate calls required during
a self energy calculation by taking advantage of the symmetries between basis integrals in an entirely
generic way

Sep 2016

*/

#include "write_tsil_ini.hpp"

//#define DEBUG

vector<int> Print_dotsil::get_duplicates()
{
  vector<int> duplicates;
  duplicates.resize(names.size());
  
  for (unsigned int j=0;j<eval_vec.size();j++)
  {
    vector<bool> tmp = V_check_vec[j];
    for (unsigned int i=0;i<names.size();i++)
    {
      if (tmp[i]){ duplicates[i] = duplicates[i] + 1;}
    }
  }
  return duplicates;
}

void Print_dotsil::print_total(std::vector<int> total)
{
  cout << " ------------------ " << endl;
  for (unsigned int j = 0; j<total.size();j++)
  {
    Bases base_tmp = base_map[names[j]];
    if (base_tmp.type == "A"){cout << total[j] << "  ";}
    else if (base_tmp.type == "B"){cout << total[j] << "   ";}
    else if (base_tmp.type == "F"){cout << total[j] << "      ";}
    else if (base_tmp.type == "V"){cout << total[j] << "     ";}
    else if (base_tmp.type == "J" || base_tmp.type == "K"|| base_tmp.type == "T"){cout << total[j] << "    ";}
    else {cout << total[j] << " - ";}
    
  }
  cout << endl;
}



void Print_dotsil::print_to_file(ofstream &myfile)
{
  cout << "basis integral sorting started" << endl;

  get_masses();
  for (unsigned int i = 0; i<names.size();i++)
  {
    get_poss_eval(base_map[names[i]]);
  }
  
  
  cout << "number of TSIL evaluate calls required before optimisation = " << eval_count << endl;
  V_check_vec.resize(eval_vec.size());
  // remove all integrals for which duplicates is less than 1
  vector<int>  sum;
  sum.resize(names.size());
  int sum_min=0;
  cout << "\n";
  int status=0;
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    vector<bool> check_vec = eval_vec[i].get_check_vec(names, base_map,masses);
    V_check_vec[i]=check_vec;
    status=(float(i)/eval_vec.size())*100;
    cout<< "\r" << "sorting integrals . . . " << status << "% complete ";
    std::cout << std::flush;
  }
  status=100;
  cout<< "\r" << "sorting integrals . . . " << status << "% complete ";
  cout << "\n";
 
 
 // remove objects that have check_vec = 0 for all values // is this possible?
  /*prestart:;
  for (unsigned int j=0;j<eval_vec.size();j++)
  {
    vector<bool> tmp = V_check_vec[j];
    bool zero = true;
    for (unsigned int i=0;i<names.size();i++)
    {
      if (tmp[i]){ zero=false;}
    }
    if (zero)
    {
      eval_vec.erase(eval_vec.begin()+j);
      V_check_vec.erase(V_check_vec.begin()+j);
      goto prestart;
    }
  }
  */
 
  
  start:;
  
  vector<int> duplicates = get_duplicates();
  for (unsigned int j=0;j<eval_vec.size();j++)
  {
    vector<bool> tmp = V_check_vec[j];
    sum_min = 0;
    for (unsigned int i=0;i<names.size();i++)
    {
      if (tmp[i])
      {
        sum[i] = duplicates[i]; sum_min = sum[i];
      }
      else { sum[i] = 0;}
    }
    
    for (unsigned int i=0;i<names.size();i++)
    {
      if (sum[i]<sum_min && tmp[i])
      {
        sum_min=sum[i];
      }
    }
    
    
    if (sum_min > 1)
    {
      // can safely remove this eval_obj
      eval_vec.erase(eval_vec.begin()+j);
      V_check_vec.erase(V_check_vec.begin()+j);
      goto start;
    }
  }
  
  
  cout << "number of TSIL evaluate calls required after optimisation = " << eval_vec.size() << endl;

  // print out matrix of check vecs that are remaining
  
  eval_obj eo_tmp_test = eval_vec[0];
  vector<bool> check_vec_tmp_test = eo_tmp_test.get_check_vec(names, base_map,masses);
  
  #ifdef DEBUG
  cout << "Size of names = " << names.size() << endl;
  cout << " -----------  check vectors are -----------  " << endl;
  for (unsigned int j = 0; j<names.size();j++)
  {
    cout << names[j] << " ";
  }
  cout << endl;
  #endif
  
  vector<int> total;
  total.resize(names.size());
  
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    eval_obj eo_tmp = eval_vec[i];
    vector<bool> check_vec_tmp = eo_tmp.get_check_vec(names, base_map,masses);
    
    for (unsigned int j = 0; j<check_vec_tmp.size();j++)
    {
      Bases base_tmp = base_map[names[j]];
      #ifdef DEBUG
      if (base_tmp.type == "A"){cout << check_vec_tmp[j] << "  ";}
      else if (base_tmp.type == "B"){cout << check_vec_tmp[j] << "   ";}
      else if (base_tmp.type == "F"){cout << check_vec_tmp[j] << "      ";}
      else if (base_tmp.type == "V"){cout << check_vec_tmp[j] << "     ";}
      else if (base_tmp.type == "J" || base_tmp.type == "K"|| base_tmp.type == "T"){cout << check_vec_tmp[j] << "    ";}
      else {cout << check_vec_tmp[j] << " - ";}
      #endif
      if (check_vec_tmp[j]){ total[j] = 1; };
      
    }
    #ifdef DEBUG
    cout << endl;
    #endif
  }
  // print total line at bottom
  #ifdef DEBUG
  print_total(total);
  #endif
  
  
  
  
  
  // print out TSIL evaluate statements
  
  
  // do the lower dimensional functions first and set total[i]=2.
  for (unsigned int i = 0;i< names.size();i++)
  {
    string type = base_map[names[i]].type;
    if (type == "A"||type=="B"||type=="K")
    {
      print_doTSIL(myfile,base_map[names[i]]);
      total[i]=2;
      if (type == "B" && base_map[names[i]].e1!=base_map[names[i]].e2)
      {
        Bases base_temp_B;
        base_temp_B.type = "B";
        base_temp_B.e1 = base_map[names[i]].e2;
        base_temp_B.e2 = base_map[names[i]].e1;
        vector<string> mass_in = {base_map[names[i]].e1,base_map[names[i]].e2};
        vector<string> id_in ={base_map[names[i]].id1,base_map[names[i]].id2};
        string name_B = get_short_name(base_temp_B,mass_in,id_in);
        myfile << name_B << " = " << base_map[names[i]].short_name << ";" << endl;
      }
      
      
    }
  }
  #ifdef DEBUG
  print_total(total);
  #endif
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    eval_obj eo_tmp = eval_vec[i];
    print_eval_obj(myfile,eo_tmp,total);
  }
  
  #ifdef DEBUG
  print_total(total);
  #endif
  
    // do the lower dimensional functions first and set total[i]=2.
  for (unsigned int i = 0;i< names.size();i++)
  {
    string type = base_map[names[i]].type;
    if (type=="J" && total[i]!=2)
    {
      print_doTSIL(myfile,base_map[names[i]]);
      
    }
  }
  
  #ifdef DEBUG
  print_total(total);
  #endif
  

}

string Print_dotsil::coeff(string type)
{
  if (type=="V"){return " - ";}
  if (type=="T"){return " - ";}
  if (type=="B"){return "i * ";}
  else {return "";}
}



void Print_dotsil::print_eval_obj(ofstream &myfile,eval_obj &eo, vector<int> &total)
{
  myfile << "TSIL_SetParameters (&bar," << eo.x << "2, " << eo.y << "2, " << eo.z << "2 , " << eo.u << "2 , " << eo.v  << "2, Q2);" << endl;
  myfile << "TSIL_Evaluate (&bar, s);" << endl;
  
  vector<bool> check_vec = eo.get_check_vec(names, base_map,masses);
  
  for (unsigned int i=0;i<names.size();i++)
  {
    if (check_vec[i] && total[i]==1)
    {
      myfile << base_map[names[i]].short_name << "="<< coeff(base_map[names[i]].type) << " TSIL_GetFunction (&bar,\""<< eo.eval_string[i] <<"\");"<< endl;
      total [i] = 2;
    }
  }
  
}



void Print_dotsil::add_eval_obj( string x, string y, string z, string u, string v)
{
  eval_obj eo(x,y,z,u,v);
  eval_vec.push_back(eo);
}


void Print_dotsil::get_poss_eval(Bases base)
{
  string a = base.e1;
  string b = base.e2;
  string c = base.e3;
  string d = base.e4;
  string e = base.e5;
  string o = "blank";

  // the following are all the possible objects that could evaluate this basis integral
  
  if (base.type == "V")
  {
    eval_count = eval_count +1;
    for (unsigned int i = 0; i<masses.size();i++)
    {
      o = masses[i];
      add_eval_obj( a , c , d , o , b );
      /*add_eval_obj( c , a , o , d , b );
      add_eval_obj( d , o , a , c , b );
      add_eval_obj( o , d , c , a , b );
      add_eval_obj( a , b , d , o , c );
      add_eval_obj( b , a , o , d , c );
      add_eval_obj( d , o , a , b , c );
      add_eval_obj( o , d , b , a , c );*/
    }
  }
  
  if (base.type == "F")
  {
    add_eval_obj( a, b, c, d, e);
    eval_count = eval_count +1;
  }
  
  


}


void eval_obj::swap_tsil_to_tarcer_V(string &a,string &b,string &c, string &d)
{
  string a_tmp = a, b_tmp = b, c_tmp = c, d_tmp = d;
  d = a_tmp;
  a = b_tmp;
  c = c_tmp;
  b = d_tmp;
}



void eval_obj::add_integral(string type,string tsil_id, string x, string y, string z, string u, string v)
{
  if (type == "V")
  {
   //cout << "--------------" << endl;
   //cout << "adding the integral (original) " << type << "  " << x << " "<< y << " "  << z << " " << u << " " << v << endl;
   swap_tsil_to_tarcer_V(x,y,z,u);
   //cout << "adding the integral (TARCER) " << type << "  " << x << " "<< y << " "  << z << " " << u << " " << v << endl;
  //cout << "--------------" << endl;
  }

  Bases base(type, x, y, z, u ,v);
  
  base.short_name = tsil_id;
  integrals.push_back(base);
}


std::vector<Bases> eval_obj::get_integrals(std::vector<string> masses_input)
{
  masses.clear();
  masses = masses_input;
  //reset integrals to a null vector
  integrals.clear();
  
  add_integral("F","M", x, y, z, u ,v);  // Master integral -- M in TSIL
  
  add_integral("V","Uzxyv", z, x, y, v); // using TARCER definition of V -- U in TSIL
  // this is a TSIL type not a FA type, need to convert //  in TARCER notation this is V x  v  y  z
  
  add_integral("V","Uuyxv", u, y ,x, v);
  /*add_integral("V","Uxzuv", x, z ,u, v);
  add_integral("V","Uyuzv", y, u, z, v);
  // switch last two arguments due to symmetry
  add_integral("V","Uzxvy", z, x, v, y);
  add_integral("V","Uuyvx", u, y ,v, x);
  add_integral("V","Uxzvu", x, z ,v, u);
  add_integral("V","Uyuvz", y, u, v, z);*/
  
  // T integrals
  
  add_integral("T","Tvyz", v, y, z );
  add_integral("T","Tuxv", u, x, v );
  add_integral("T","Tyzv", y, z, v );
  add_integral("T","Txuv", x, u, v );
  add_integral("T","Tzyv", z, y, v );
  add_integral("T","Tvxu", v, x, u );
  // permutations of T integrals are last two arguments swapped
  add_integral("T","Tvzy", v, z, y );
  add_integral("T","Tuvx", u, v, x );
  add_integral("T","Tyvz", y, v, z );
  add_integral("T","Txvu", x, v, u );
  add_integral("T","Tzvy", z, v, y );
  add_integral("T","Tvux", v, u, x );
  
  add_integral("B","Bxz",x,z);
  add_integral("B","Bzx",z,x);
  add_integral("B","Byu",y,u);
  add_integral("B","Buy",u,y);
  
  // "S" integrals are what we refer to as "J" integrals
  
  add_integral("J","Svyz",v,y,z);
  add_integral("J","Suxv",u,x,v);
  // symmetries give:
  add_integral("J","Sxuv",x, u, v);
  add_integral("J","Suvx",u, v, x);
  add_integral("J","Svux",v, u, x);
  add_integral("J","Sxvu",x, v, u);
  add_integral("J","Svxu",v, x, u);
  
  add_integral("J","Svzy",v,z,y);
  add_integral("J","Szyv",z,y,v);
  add_integral("J","Szvy",z,v,y);
  add_integral("J","Syvz",y,v,z);
  add_integral("J","Syzv",y,z,v);
  
  return integrals;
  
}


std::vector<bool> eval_obj::get_check_vec(vector<string> names, std::map<std::string, Bases> base_map,std::vector<string> masses_input)
{
  if (check_vec_evaluated){return check_vec;}
  else
  {
    vector<Bases> basis_objects = get_integrals(masses_input);
    check_vec.clear();
    check_vec.resize(names.size());
    location.resize(names.size());
    eval_string.resize(names.size());
    
    for (unsigned int i = 0; i<basis_objects.size();i++)
    {
      Bases b_tmp = basis_objects[i];
      
      for (unsigned int j = 0; j<names.size();j++)
      {
        Bases basis = base_map[names[j]];
        if ( (b_tmp.type == basis.type) && (b_tmp.e1 == basis.e1) && (b_tmp.e2 == basis.e2) && (b_tmp.e3 == basis.e3) && (b_tmp.e4 == basis.e4)  && (b_tmp.e5 == basis.e5))
        {
          check_vec[j] = true;
          location[j] = i;
          eval_string[j] = b_tmp.short_name;
        }
      }
    }
    check_vec_evaluated = true;
    return check_vec;
  }
}



