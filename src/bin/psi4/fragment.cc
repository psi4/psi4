#include "fragment.h"
#include "element_to_Z.h"

namespace psi { namespace opt09 {

extern "C" { extern FILE *outfile; }

// deep copy constructor
Fragment::Fragment(const Fragment & frag)
{
  int i,j;
  natoms = frag.natoms;

  Z = new double [natoms];
  for (i=0; i<natoms; ++i) Z[i] = frag.Z[i];

  masses = new double [natoms];
  for (i=0; i<natoms; ++i) masses[i] = frag.masses[i];

  atom_label = new string [natoms];
  for (i=0; i<natoms; ++i)
    atom_label[i] = frag.atom_label[i];

  if (frag.geom != NULL) {
    geom = block_matrix(natoms,3);
    for (i=0; i<natoms; ++i)
      for (j=0; j<3; ++j) 
        geom[i][j] = frag.geom[i][j];
  }
}

// deep assignment operator
Fragment & Fragment::operator=(const Fragment & frag)
{
  if (this == &frag) //object assigned to itself
    return *this;

  // delete lhs
  delete [] Z;
  delete [] masses;
  delete [] atom_label;
  if (geom != NULL) { free_block(geom); geom = NULL; }

  // copy rhs
  int i,j;
  natoms = frag.natoms;

  Z = new double [natoms];
  for (i=0; i<natoms; ++i) Z[i] = frag.Z[i];

  masses = new double [natoms];
  for (i=0; i<natoms; ++i) masses[i] = frag.masses[i];

  atom_label = new string [natoms];
  for (i=0; i<natoms; ++i) 
    atom_label[i] = frag.atom_label[i];

  if (frag.geom != NULL) {
    geom = block_matrix(natoms,3);
    for (i=0; i<natoms; ++i)
      for (j=0; j<3; ++j)
        geom[i][j] = frag.geom[i][j];
  }
  return *this;
}

// does not support non-"simple" geometry format:
// geometry = ( (atom 1) (atom 2) ... )
void Fragment::read_cartesian_from_input(string geom_label,double conv_factor)
{
  int i, xyz, errcod, num_elem;
  double tmp;
  char *c_geom_label, *c_atom_label;
  string error_message, string_atom_label;

  c_geom_label = const_cast<char *>(geom_label.c_str());

  // Determine number of atoms in fragment: natoms
  num_elem = 0;
  ip_count(c_geom_label, &num_elem, 0);
  if (!natoms) {
    error_message = "Fragment " + geom_label + " has no atoms.";
    throw(error_message);
  }
  if (num_elem % 4) {
    error_message = "Problem with number of entries in " + geom_label;
    throw(error_message);
  }
  natoms = num_elem/4;

  // allocate memory for Z, geom
  Z = new double [natoms];
  masses = new double [natoms];
  atom_label = new string [natoms]; 
  geom = block_matrix(natoms,3);

  for (i=0; i<natoms; ++i) {
    errcod = ip_string(c_geom_label,&c_atom_label,1,4*i);
    if (errcod != IPE_OK) {
      error_message = "Problem reading the " + geom_label + " array.";
      throw(error_message);
    }
    Element_to_Z elem_map;
    Z[i] = elem_map[string_atom_label = c_atom_label]; // convert to string
    atom_label[i] = string_atom_label;
    free(c_atom_label);

    for (xyz=0; xyz<3; ++xyz) {
      tmp = 0.0;
      errcod = ip_data(c_geom_label,"%lf", &tmp,1,4*i+xyz+1);
      if (errcod != IPE_OK) {
        error_message = "Problem reading the " + geom_label + " array.";
        throw(error_message);
      }
      geom[i][xyz] = tmp;
    }
  }
  for (i=0; i<natoms; ++i)
    for (xyz=0; xyz<3; ++xyz)
      geom[i][xyz] *= conv_factor;

  return;
}

// read masses from input file - return 1 if successful
int Fragment::read_masses_from_input(int fragment_id) throw(bad_masses_io)
{
  int i,j,a,cnt;
  char *buf, *c_mass_label;
  double tval;
  string mass_label;
  bool found;

  ostringstream outstr;
  if (fragment_id == 0)
    outstr << "ISOTOPES";
  else
    outstr << "ISOTOPES" << fragment_id;

  mass_label = outstr.str();
  outstr.clear();
  c_mass_label = const_cast<char *>(mass_label.c_str());

  cnt = 0;
  if (ip_exist(c_mass_label,0)) {
    a = 0;
    ip_count(c_mass_label, &a, 0);
    if (a != natoms) {
      fprintf(outfile,"%s array has wrong dimension.\n", c_mass_label);
      throw bad_masses_io(c_mass_label,0);
    }
    for (i=0;i<natoms;++i) {
      found = false;
      ip_string(c_mass_label,&buf,1,i);
      for (j=0;j<LAST_MASS_INDEX;j++) {
        if (!strcmp(buf, mass_labels[j])) {
          masses[cnt++] = atomic_masses[j];
          found = 1;
          break;
        }
      }
      if (!found) {
        fprintf(outfile, "Isotope label %s is unidentifiable.\n", buf);
        throw bad_masses_io(c_mass_label,i);
      }
    }
    return 1;
  }


  if (fragment_id == 0) 
    outstr << "MASSES";
  else
    outstr << "MASSES" << fragment_id;

  mass_label = outstr.str();
  outstr.clear();
  c_mass_label = const_cast<char *>(mass_label.c_str());

  if (ip_exist(c_mass_label,0)) {
    a = 0;
    ip_count(c_mass_label,&a,0);
    if (a != natoms) {
      fprintf(outfile,"%s array has wrong dimension\n", c_mass_label);
      throw bad_masses_io(c_mass_label,0);
    }
    else {
      for(i=0;i<natoms;++i) {
        ip_data(c_mass_label,"%lf",&tval,1,i);
        if (tval < 0.0) {
          fprintf(outfile,"Given mass value is negative!\n");
          throw bad_masses_io(c_mass_label,i);
        }
        masses[cnt++] = tval;
      }
    }
    return 1;
  }

  return 0;
}

// use atomic numbers to set masses
void Fragment::read_default_masses(void)
{
  int cnt=0;
  for(int i=0; i<natoms; ++i)
    masses[cnt++] = an2masses[(int) Z[i]];
}

}}

