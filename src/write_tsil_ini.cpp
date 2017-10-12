/*
 Mass Builder
 
 James McKay
 Sep 2016
 
 --- write_tsil_ini.cpp ---
 
 The functions defined in this file add the ability to significantly reduce the number TSIL_Evaluate calls required during
 a self energy calculation by taking advantage of the symmetries between basis integrals in an entirely
 generic way
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
  
  cout << "number of full TSIL evaluate calls required before optimisation = " << eval_count << endl;
  V_check_vec.resize(eval_vec.size());
  // remove all integrals for which duplicates is less than 1
  vector<int>  sum;
  sum.resize(names.size());
  int sum_min=0;
  cout << "\n";
  int status=0;
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    vector<bool> check_vec = eval_vec[i].get_check_vec(names, base_map);
    V_check_vec[i]=check_vec;
    
    
    status=(float(i)/eval_vec.size())*100;
    
    std::string faces[7] = {":-C",":-(",":-|",":-)",":-D",":-O"};
    
    int fc = floor((status * 7.0)/100);
    
	  cout<< "\r" << "sorting integrals . . . " << status << "% complete " << faces[fc+1];
	  std::cout << std::flush;
    
  }
  status=100;
  cout<< "\r" << "sorting integrals . . . " << status << "% complete :-O";
  cout << "\n";
  
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
  
  
  cout << "number of full TSIL evaluate calls required after optimisation = " << eval_vec.size() << endl;
  
  // print out matrix of check vecs that are remaining
  
  eval_obj eo_tmp_test = eval_vec[0];
  vector<bool> check_vec_tmp_test = eo_tmp_test.get_check_vec(names, base_map);
  
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
    vector<bool> check_vec_tmp = eo_tmp.get_check_vec(names, base_map);
    
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
  
  // option to not evaluate the following two-loop integrals
  myfile << "  if (data.do_tsil_all)" << endl;
  myfile << "  {" << endl;
  
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    eval_obj eo_tmp = eval_vec[i];
    print_eval_obj(myfile,eo_tmp,total);
  }
  for (unsigned int i = 0;i< names.size();i++)
  {
    string type = base_map[names[i]].type;
    if (type=="J" && total[i]!=2)
    {
      print_doTSIL(myfile,base_map[names[i]]);
      
    }
  }
  myfile << "  }" << endl;
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
  myfile << "    TSIL_SetParameters (&bar," << eo.x << "2, " << eo.y << "2, " << eo.z << "2 , " << eo.u << "2 , " << eo.v  << "2, Q2);" << endl;
  myfile << "    TSIL_Evaluate (&bar, s);" << endl;
  
  vector<bool> check_vec = eo.get_check_vec(names, base_map);
  
  for (unsigned int i=0;i<names.size();i++)
  {
    if (check_vec[i] && total[i]==1)
    {
      myfile << "    " << base_map[names[i]].short_name << "="<< coeff(base_map[names[i]].type) << " TSIL_GetFunction (&bar,\""<< eo.eval_string[i] <<"\");"<< endl;
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
      add_eval_obj( d , c , b , o , a );
      add_eval_obj( d , c , b , o , a );
      add_eval_obj( c , d , o , b , a );
      add_eval_obj( b , o , d , c , a );
      add_eval_obj( o , b , c , d , a );
      add_eval_obj( d , a , b , o , c );
      add_eval_obj( a , d , o , b , c );
      add_eval_obj( b , o , d , a , c );
      add_eval_obj( o , b , a , d , c );
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
  a = d_tmp;
  b = a_tmp;
  c = c_tmp;
  d = b_tmp;
}



void eval_obj::add_integral(string type,string tsil_id, string a, string b, string c, string d, string e)
{
  if (type == "V")
  {
    swap_tsil_to_tarcer_V(a,b,c,d);
  }
  
  Bases base(type, a, b, c, d, e);
  
  base.short_name = tsil_id;
  integrals.push_back(base);
}


std::vector<Bases> eval_obj::get_integrals()
{
  integrals.clear();
  
  // F (TARCER) / M (TSIL) integral
  add_integral("F","M", x, y, z, u ,v);
  
  // V (TARCER) / U (TSIL) integrals
  add_integral("V","Uzxyv", z, x, y, v);
  add_integral("V","Uuyxv", u, y ,x, v);
  add_integral("V","Uxzuv", x, z ,u, v);
  add_integral("V","Uyuzv", y, u, z, v);
  // switch last two arguments due to symmetry
  add_integral("V","Uzxvy", z, x, v, y);
  add_integral("V","Uuyvx", u, y ,v, x);
  add_integral("V","Uxzvu", x, z ,v, u);
  add_integral("V","Uyuvz", y, u, v, z);
  
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
  
  // Tbar integrals
  
  add_integral("TBAR","TBARvyz", v, y, z );
  add_integral("TBAR","TBARuxv", u, x, v );
  add_integral("TBAR","TBARyzv", y, z, v );
  add_integral("TBAR","TBARxuv", x, u, v );
  add_integral("TBAR","TBARzyv", z, y, v );
  add_integral("TBAR","TBARvxu", v, x, u );
  // permutations of T integrals are last two arguments swapped
  add_integral("TBAR","TBARvzy", v, z, y );
  add_integral("TBAR","TBARuvx", u, v, x );
  add_integral("TBAR","TBARyvz", y, v, z );
  add_integral("TBAR","TBARxvu", x, v, u );
  add_integral("TBAR","TBARzvy", z, v, y );
  add_integral("TBAR","TBARvux", v, u, x );  
  
  add_integral("B","Bxz",x,z);
  add_integral("B","Bzx",z,x);
  add_integral("B","Byu",y,u);
  add_integral("B","Buy",u,y);
  
  // TSIL "S" integrals are what we refer to as "J" integrals from TARCER
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


std::vector<bool> eval_obj::get_check_vec(vector<string> names, std::map<std::string, Bases> base_map)
{
  if (check_vec_evaluated){return check_vec;}
  else
  {
    vector<Bases> basis_objects = get_integrals();
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



