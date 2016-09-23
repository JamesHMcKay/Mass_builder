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

double iter_factorial(int n)
{
    double result = 1;
    for(int i = 1; i <= n; i++)
    {
        result = result * i;
      
    }
    return result;
}

void Print_dotsil::make_sets()
{
    vector<eval_obj> subset;
    vector< vector<eval_obj> > all_sets;
    int n,len_in;// n is the number of objects in each set, len is the number of eval_objects
    len_in = eval_vec.size()-1;
    cout << "len_in = " << len_in << endl;
    cin>>n;
    std::vector<int> a, len,i;
    a.resize(n);
    len.resize(n);
    //cin>>len_in;
    int j = 999;
    for(int i = 0 ; i < n ; i++)
    {
        len[i]=len_in;
        a[i]=0;
    }
  
    while(1)
    {
        subset = {};
        vector<int> b = a;
        sort(b.begin(),b.end());
        b.erase( unique( b.begin(), b.end() ), b.end() );
        if (b.size()==n)
        {
          for(unsigned int j = 0 ; j<b.size();j++)
          {
            subset.push_back( eval_vec[b[j]] );
            cout<<b[j]<<" ";
          }
          all_sets.push_back(subset);
          cout<<endl;
        }
        for(j = n-1 ; j>=0 ; j--)
        {
            if(++a[j]<=len[j])
            {
                for (int k=j;k<n;k++)
                {
                  a[k] = a[j];
                }
                break;

            }
            else
            {
                a[j]=0;
             }
        }
        if(j<0)
        {
            break;
        }
    }
  cout << "length of unique elements = " << all_sets.size() << endl;
  cout << "from set of " << len_in + 1 << " elements we choose " << n << endl;
  cout << "this should be equal to = " << iter_factorial(len_in + 1) / ( iter_factorial(n) * iter_factorial(len_in + 1-n)) << endl;
}


vector<int> Print_dotsil::get_duplicates()
{
  cout  << "evaluating duplicates " << endl;
  vector<int> duplicates;
  duplicates.resize(names.size());
  for (unsigned int i=0;i<names.size();i++)
  {
    duplicates[i]=0;
    for (unsigned int j=0;j<eval_vec.size();j++)
    {
      vector<bool> tmp = V_check_vec[j];
      if (tmp[i]){ duplicates[i] = duplicates[i] + 1;}
      
    }
   
  }
  return duplicates;
}


void Print_dotsil::sort_integrals()
{

  // testing

  /*Bases basis("V","Mb","Mb","Mb","Mb");
  base_map["int1"] = basis;
  Bases basis2("V","Mb","Mb","Mc","Mb");
  base_map["int2"] = basis2;
  Bases basis3("V","Mb","Mb","Mb","Mg");
  base_map["int3"] = basis3;
  Bases basis4("F","Mb","Mb","Mb","Mg","Mb");
  base_map["int4"] = basis4;
  names = {"int1","int2","int3","int4"};
  */
  
  get_masses();
  for (unsigned int i = 0; i<names.size();i++)
  {
    get_poss_eval(base_map[names[i]]);
  }
  
  
  cout << "number of eval objects = " << eval_vec.size() << endl;
  V_check_vec.resize(eval_vec.size());
  // remove all integrals for which duplicates is less than 1
  vector<int>  sum;
  sum.resize(names.size());
  int sum_min=0;
  
  start:;
  
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    eval_obj eo_tmp = eval_vec[i];
    vector<bool> check_vec = eo_tmp.get_check_vec(names, base_map,masses);
    V_check_vec[i]=check_vec;
  }
  vector<int> duplicates = get_duplicates();
  for (unsigned int j=0;j<eval_vec.size();j++)
  {
    vector<bool> tmp = V_check_vec[j];
    sum_min = 0;
    for (unsigned int i=0;i<names.size();i++)
    {
    if (tmp[i]){ sum[i] = duplicates[i];sum_min = sum[i];}
    else { sum[i] = 0;}
    }
    // find minimum element
    
    for (unsigned int i=0;i<names.size();i++)
    {
    if (sum[i]<sum_min && tmp[i]){sum_min=sum[i];}
    }
    
    
    if (sum_min > 1)
    {
    // can safely remove this eval_obj
    eval_vec.erase(eval_vec.begin()+j);
    V_check_vec.erase(V_check_vec.begin()+j);
    goto start;
    }
  }
  
  // of the remaining sets keep
  
  
  
  
  
  // print out remaing eval statements
  
  for (unsigned int i = 0; i<eval_vec.size();i++)
  {
    eval_obj eo_tmp = eval_vec[i];
    cout << "eval object = " << eo_tmp.x << " " << eo_tmp.y << " " << eo_tmp.z << " "<< eo_tmp.u << " " << eo_tmp.v << endl;
    std::vector<Bases> integrals_tmp = eo_tmp.get_integrals(masses);
    vector<bool> check_vec_tmp = eo_tmp.get_check_vec(names, base_map,masses);
    
    for (unsigned int j = 0; j<check_vec_tmp.size();j++)
    {
      if (check_vec_tmp[j])
      {
        Bases base_tmp = integrals_tmp[eo_tmp.location[j]];
        cout << "print integral " << base_tmp.type << "( " << base_tmp.e1 << " , " << base_tmp.e2 << " , "  << base_tmp.e3 << " , " << base_tmp.e4 << " , " << base_tmp.e5 << " )" << endl;
      }
    
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
    for (unsigned int i = 0; i<masses.size();i++)
    {
      o = masses[i];
      add_eval_obj( o, a, d, b, c);
      add_eval_obj( a, o, b, c, d);
      add_eval_obj( a, o, b, c, d);
      add_eval_obj( o, a, c, b, d);
      add_eval_obj( b, d, a, o, c);
      add_eval_obj( d, b, o, a, c);
      add_eval_obj( b, c, a, o, d);
      add_eval_obj( a, o, b, d, c);
    }
  }
  
  if (base.type == "F")
  {
    add_eval_obj( a, b, c, d, e);
  }


}


void eval_obj::add_integral(string type, string x, string y, string z, string u, string v)
{
  string xi=x, yi=y, zi=z, ui=u, vi=v;
  for (unsigned int i = 0; i<masses.size();i++)
  {
    if (x == "blank"){ xi = masses[i];}
    if (y == "blank"){ yi = masses[i];}
    if (z == "blank"){ zi = masses[i];}
    if (u == "blank"){ ui = masses[i];}
    if (v == "blank"){ vi = masses[i];}
    Bases base(type, xi, yi, zi, ui ,vi);
    integrals.push_back(base);
  }
}



std::vector<Bases> eval_obj::get_integrals(std::vector<string> masses_input)
{
  masses.clear();
  masses = masses_input;
  //reset integrals to a null vector
  integrals.clear();
  
  add_integral("F", x, y, z, u ,v);  // Master integral -- M in TSIL
  add_integral("V", z, x, y, v, ""); // using TARCER definition of V -- U in TSIL
  add_integral("V", u, y ,x, v, "");
  add_integral("V", x, z ,u, v, "");
  add_integral("V", y, u, z, v, "");
  // switch last two arguments due to symmetry
  add_integral("V", z, x, v, y);
  add_integral("V", u, y ,v, x);
  add_integral("V", x, z ,v, u);
  add_integral("V", y, u, v, z);
  
  // T integrals
  
  add_integral("T", v, y, z );
  add_integral("T", u, x, v );
  add_integral("T", y, z, v );
  add_integral("T", x, u, v );
  add_integral("T", z, y, v );
  add_integral("T", v, x, u );
  // permutations of T integrals are last two arguments swapped
  add_integral("T", v, z, y );
  add_integral("T", u, v, x );
  add_integral("T", y, v, z );
  add_integral("T", x, v, u );
  add_integral("T", z, v, y );
  add_integral("T", v, u, x );
  
  add_integral("B",x,z);
  add_integral("B",z,x);
  add_integral("B",y,u);
  add_integral("B",u,y);
  
  // "S" integrals are what we refer to as "J" integrals
  
  add_integral("J",v,y,z);
  add_integral("J",u,x,v);
  // symmetries give:
  add_integral("J", x,u ,v );
  add_integral("J", u, v, x);
  add_integral("J", v, u, x);
  add_integral("J", x, v, u);
  add_integral("J", v, x, u);
  
  add_integral("J",v,z,y);
  add_integral("J",z,y,v);
  add_integral("J",z,v,y);
  add_integral("J",y,v,z);
  add_integral("J",y,z,v);
  
  return integrals;
  
}


std::vector<bool> eval_obj::get_check_vec(vector<string> names, std::map<std::string, Bases> base_map,std::vector<string> masses_input)
{
  vector<Bases> basis_objects = get_integrals(masses_input);
  check_vec.clear();
  check_vec.resize(names.size());
  location.resize(names.size());
  
  
  
  for (unsigned int j = 0; j<names.size();j++)
  {
    Bases basis = base_map[names[j]];
    for (unsigned int i = 0; i<basis_objects.size();i++)
    {
    Bases b_tmp = basis_objects[i];
    if ( (b_tmp.type == basis.type) && (b_tmp.e1 == basis.e1) && (b_tmp.e2 == basis.e2) && (b_tmp.e3 == basis.e3) && (b_tmp.e4 == basis.e4))
    {
      check_vec[j] = true;
      location[j] = i;
    }
    }
  }



return check_vec;


}



