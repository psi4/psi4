#include <Fragment.h>
#include <Element_to_Z.h>

namespace psi {

using std::string;

// deep copy constructor
Fragment::Fragment(const Fragment & frag)
{
  int i,j;
  num_atoms = frag.num_atoms;
  Z = init_array(num_atoms);
  geom = block_matrix(num_atoms,3);
  atom_label = new string [num_atoms];
  for (i=0; i<num_atoms; ++i) {
    Z[i] = frag.Z[i];
    atom_label[i] = frag.atom_label[i];
    for (j=0; j<3; ++j) 
      geom[i][j] = frag.geom[i][j];
  }
  if (frag.masses != NULL) {
    masses = init_array(num_atoms);
    for (i=0; i<num_atoms; ++i)
      masses[i] = frag.masses[i];
  }
}

// deep assignment operator
Fragment & Fragment::operator=(const Fragment & frag)
{
  if (this == &frag) //object assigned to itself
    return *this;

  delete [] Z;   // delete old fragment
  delete [] atom_label;
  if (geom != NULL) free_block(geom);

  int i,j; //deep copy fragment
  num_atoms = frag.num_atoms;
  geom = block_matrix(num_atoms,3);
  atom_label = new string [num_atoms];
  for (i=0; i<num_atoms; ++i) {
    Z[i] = frag.Z[i];
    atom_label[i] = frag.atom_label[i];
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

  //if (Z != NULL) delete [] Z;  // clear current contents, if any
  //if (atom_label != NULL) delete [] atom_label;
  //if (geom != NULL) free_block(geom);

  c_geom_label = const_cast<char *>(geom_label.c_str());

  // Determine number of atoms in fragment
  num_elem = 0;
  ip_count(c_geom_label, &num_elem, 0);
  if (!num_atoms) {
    error_message = "Fragment " + geom_label + " has no atoms.";
    throw(error_message);
  }
  if (num_elem % 4) {
    error_message = "Problem with number of entries in " + geom_label;
    throw(error_message);
  }
  num_atoms = num_elem/4;

  // allocate memory for Z, geom
  Z = init_array(num_atoms);
  geom = block_matrix(num_atoms,3);
  atom_label = new string [num_atoms]; 

  for (i=0; i<num_atoms; ++i) {
    // read symbol and determine atomic number
    errcod = ip_string(c_geom_label,&c_atom_label,1,4*i);
    if (errcod != IPE_OK) {
      error_message = "Problem reading the " + geom_label + " array.";
      throw(error_message);
    }
    Element_to_Z elem_map;
    Z[i] = elem_map[string_atom_label = c_atom_label];
    atom_label[i] = string_atom_label;
    free(c_atom_label);

    // read geometry
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
  for (i=0; i<num_atoms; ++i)
    for (xyz=0; xyz<3; ++xyz)
      geom[i][xyz] *= conv_factor;

  return;
}

// read masses from input file - return 1 if successful
int Fragment::read_masses_from_input(string geom_label) throw(bad_masses_io)
{
  int i,j,a,cnt,natom;
  char *buf;
  double tval;

  char *c_geom_label;
  c_geom_label = const_cast<char *>(geom_label.c_str());
  natom = get_natom();

  if (masses == NULL) masses = new double[natom];

  cnt = 0;
  if (ip_exist("ISOTOPES",0)) {
    a = 0;
    ip_count("ISOTOPES", &a, 0);
    if (a != natom) {
      fprintf(outfile,"ISOTOPES array has wrong dimension.\n");
      throw bad_masses_io("ISOTOPES",0);
    }
    for (i=0;i<natom;++i) {
      ip_data("ISOTOPES","%s", &buf,1,i);
      for (j=0;j<LAST_MASS_INDEX;j++) {
        if (!strcmp(buf, mass_labels[j])) {
          masses[cnt++] = atomic_masses[j];
          break;
        }
      }
      fprintf(outfile, "Isotope label %s is unidentifiable.\n", buf);
      throw bad_masses_io("ISOTOPES",i);
    }
    return 1;
  }
  else if (ip_exist("MASSES",0)) {
    a = 0;
    ip_count("MASSES",&a,0);
    if (a != natom) {
      fprintf(outfile,"MASSES array has wrong dimension\n");
      throw bad_masses_io("MASSES",0);
    }
    else {
      for(i=0;i<natom;++i) {
        ip_data("MASSES","%lf",&tval,1,i);
        if (tval < 0.0) {
          fprintf(outfile,"Given mass value is negative!\n");
          throw bad_masses_io("MASSES",i);
        }
        masses[cnt++] = tval;
      }
    }
    return 1;
  }
  else
    return 0;
}

// use atomic numbers to set masses
void Fragment::load_default_masses(void)
{
  for(int i=0; i<get_natom(); ++i)
    masses[cnt++] = an2masses[(int) Z[i]];
}

} // psi

