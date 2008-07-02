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

} // psi

