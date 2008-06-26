#ifndef _psi4_src_lib_fragment_h_
#define _psi4_src_lib_fragment_h_

#include <string>

namespace psi {

using std::vector;
using std::cout;
using std::string;
using std::endl;

class Fragment {
  int num_atoms;  // number of atoms in fragment
  int *Z;         // atomic numbers
  double *masses; // atomic masses
  double **geom;   // cartesian coordinates
 public:
  Fragment() {};  // default constructor is empty
  ~Fragment() {
    delete [] Z;
    delete [] masses;
    if (geom != NULL) free_block(geom);
  }
  // read cart geom from input file and "$geom_label ="
  read_cartesian_from_input(char *geom_label);
}

Fragment::Fragment(void) {
}

Fragment::read_cartesian_from_input(char *geom_label) {
  int i, j, errcod, atomcount, f, all_atomcount;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label, **geom_label, error_message[80];
  int ival;

  int num_elem, simple_geom, entry_length;

  num_elem = 0;
  ip_count(geom_label, &num_elem, 0);

  // is geometry simple, geometry = (H 0.0 0.0 0.0 C 0.0 ...)?
  simple_geom = 1;
  entry_length = 0;
  for(i=0; i < num_elem; i++) {
    ip_count(geom_label, &entry_length, 1, i);
    if(entry_length > 1) simple_geom = 0;
  }

  if (simple_geom) {
    if (num_elem%4) {
      sprintf(error_message, "Problem with number of elements in %s.",geom_label);
      throw(error_message);
    }
    num_atoms = num_elem/4;
    if (!num_atoms)
      throw("fragment has no atoms");
  }
  else {
    num_atoms = num_elem;
  }

  geom = block_matrix(num_atoms,3);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  elemsymb_charges = init_array(num_atoms);

  for (i=0; i<num_atoms; ++i) {
    if(simple_geom)
      errcod = ip_string(geom_label,&atom_label,1,4*i);
    else
      errcod = ip_string(geom_label,&atom_label,2,i,0);

    if (errcod != IPE_OK) {
      sprintf(error_message,"Problem with the %s array.", geom_label);
      throw(error_message);
    }
    if (strcmp(atom_label,"X")) {
         atom_num(atom_label, &Z);
         elemsymb_charges[atomcount] = Z;
         element[atomcount] = elem_name[(int)Z];
         atom_dummy[all_atomcount] = 0;
         ++atomcount;
       }
       else {
         full_element[all_atomcount] = "X";
         atom_dummy[all_atomcount] = 1;
       }
       free(atom_label);

      for(j=0; j<3;j++){
        if(simple_geom)
          errcod = ip_data(geom_label,"%lf", &tmp,1,4*i+j+1);
        else
          errcod = ip_data(geom_label,"%lf", &tmp,2,i,j+1);
        if (errcod != IPE_OK) {
          sprintf(error_message,"Problem with the %s array.", geom_label);
          throw(error_message);
        }
        else
          full_geom[all_atomcount][j] = tmp*conv_factor;
      }
      ++all_atomcount;
    }
  }
  read_charges();

  return;
}

}

#endif
