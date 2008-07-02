/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <Molecular_system.h>

namespace psi { namespace input {

void read_cart(Molecular_system & molecules)
{
  int i, j, errcod, atomcount, f, all_atomcount;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label, error_message[80];
  int simple_geom, num_elem, entry_length, ival;

  num_atoms        = molecules.get_num_atoms();
  num_allatoms     = num_atoms;

  // PSI3 did NOT include dummies in these
  elemsymb_charges = molecules.get_Z(); 
  nuclear_charges = molecules.get_Z(); 
  geometry     = molecules.get_geom();

  std::string *estring = molecules.get_atom_label();

  element = (char **) malloc(num_atoms*sizeof(char *));
  for (i=0; i<num_atoms; ++i) {
    element[i] = (char *) malloc((estring[i].size()+1)*sizeof(char));
    for (j=0; j<estring[i].size(); ++j)
      element[i][j] = estring[i][j];
    element[i][estring[i].size()] = '\0';
  }

  full_element = (char **) malloc(num_atoms*sizeof(char *));
  for (i=0; i<num_atoms; ++i) {
    full_element[i] = (char *) malloc((estring[i].size()+1)*sizeof(char));
    for (j=0; j<estring[i].size(); ++j)
      full_element[i][j] = estring[i][j];
    full_element[i][estring[i].size()] = '\0';
  }

  // PSI3 did include dummies in these
  full_geom        = molecules.get_geom();
  atom_dummy   = (int *) malloc(sizeof(int)*num_allatoms);

  //std::string *atom_label;
  //atom_label =  molecules.get_atom_label();

  //read_charges();

  return;
}

}} // namespace psi::input
