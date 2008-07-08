
/*
template<typename intco_type>
void compute_vtype(vector<intco_type> & intcos)
{
  for (int i=0; i<intcos.size(); ++i)
    intcos[i].compute(geom);
}

template<typename intco_type>
void compute_s_vtype(vector<intco_type> & intcos)
{
  for (int i=0; i<intcos.size(); ++i)
    intcos[i].compute_s(geom);
}

template<typename intco_type>
void print_vtype(vector<intco_type> & intcos, int print_flag)
{
  for (int i=0; i<intcos.size(); ++i)
    intcos[i].print(print_flag);
}

template<typename intco_type>
void print_s_vtype(vector<intco_type> & intcos) {
  for (int i=0; i<intcos.size(); ++i)
    intcos[i].print_s();
}
*/

template<typename intco_type>
bool is_unique(vector<intco_type> & intcos, intco_type & newone) {
  for (int i=0; i<intcos.size(); ++i)
    if (newone == intcos.at(i))
      return false; 
  return true;
}

template<typename intco_type>
int build_B_vtype(vector<intco_type> & intcos, int row) {
  int index; // runs through atoms involved in coordinate
  int i, xyz;
  for (i=0; i<intcos.size(); ++i) {
    for (xyz=0; xyz<3; ++xyz) {
      for (index=0; index<intcos.at(i).get_natom(); ++index)
        Bsimp[row][3*intcos.at(i).get_atom(index)+xyz] =
          intcos.at(i).get_s(index,xyz);
    }
    ++row;
  }
  return row;
}
