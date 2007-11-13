/*! \file 
    \ingroup (INPUT)
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

namespace psi { namespace input {

void read_cart()
{
  int i, j, errcod;
  int atomcount;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label;
  int simple_geom, num_elem, entry_length;

  ip_count("GEOMETRY", &num_elem, 0);
  simple_geom = 1;
  entry_length = 0;
  for(i=0; i < num_elem; i++) {
    ip_count("GEOMETRY", &entry_length, 1, i);
    if(entry_length > 1) simple_geom = 0;
  }

  if(simple_geom && num_elem%4) 
    punt("Problem with number of elements in GEOMETRY.");

  if(simple_geom) {
    num_allatoms = num_elem/4;
  }
  else {
    num_allatoms = 0;
    ip_count("GEOMETRY",&num_allatoms,0);
  }

  if (num_allatoms == 0)
    punt("GEOMETRY is empty!");
  else if (num_allatoms > MAXATOM)
    punt("There are more atoms than allowed!");

  /* Figure out how many non-dummy atoms are there */
  num_atoms = 0;
  for(i=0;i<num_allatoms;i++){
    if(simple_geom) 
      errcod = ip_string("GEOMETRY", &atom_label, 1, i*4);
    else
      errcod = ip_string("GEOMETRY",&atom_label,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem reading GEOMETRY array.");
    if (strcmp(atom_label,"X"))
      ++num_atoms;
  }

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  full_geom = block_matrix(num_allatoms,3);
  geometry = (double **) malloc(num_atoms*sizeof(double *));
  atom_dummy = (int *) malloc(sizeof(int)*num_allatoms);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  full_element = (char **) malloc(sizeof(char *)*num_allatoms);
  elemsymb_charges = init_array(num_atoms);

  atomcount = 0;
  for(i=0;i<num_allatoms;i++){
    if(simple_geom) 
      errcod = ip_string("GEOMETRY",&atom_label,1,4*i);
    else
      errcod = ip_string("GEOMETRY",&atom_label,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem with the GEOMETRY array.");
    if (strcmp(atom_label,"X")) {
       atom_num(atom_label, &Z);
       free(atom_label);
       elemsymb_charges[atomcount] = Z;
       element[atomcount] = elem_name[(int)Z];
       full_element[i] = elem_name[(int)Z];
       geometry[atomcount] = full_geom[i];
       atom_dummy[i] = 0;
       ++atomcount;
     }
     else {
       full_element[i] = "X";
       free(atom_label);
       atom_dummy[i] = 1;
     }

    for(j=0; j<3;j++){
      if(simple_geom) 
        errcod = ip_data("GEOMETRY","%lf", &tmp,1,4*i+j+1);
      else
        errcod = ip_data("GEOMETRY","%lf", &tmp,2,i,j+1);
      if (errcod != IPE_OK)
	punt("Problem with the GEOMETRY array.");
      else
	full_geom[i][j] = tmp*conv_factor;
    }
  }

  read_charges();

  return;
}

}} // namespace psi::input
