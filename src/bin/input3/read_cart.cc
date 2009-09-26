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
#include <libchkpt/chkpt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void read_cart()
{
  int i, j, errcod, atomcount, f, all_atomcount;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label, **geom_label, error_message[80];
  int simple_geom, num_elem, entry_length, ival;

  frag_num_atoms = (int *) malloc(nfragments*sizeof(int));
  frag_num_allatoms = (int *) malloc(nfragments*sizeof(int));
  frag_atom = (int *) malloc(nfragments*sizeof(int));
  frag_allatom = (int *) malloc(nfragments*sizeof(int));
  nref_per_fragment = (int *) malloc(nfragments*sizeof(int));
  ref_pts_lc = (double ***) malloc(nfragments*sizeof(double **));
  num_atoms = 0;
  num_allatoms = 0;

  geom_label = (char **) malloc(nfragments*sizeof(char *));

  for (f=0; f<nfragments; ++f) {

    geom_label[f] = (char *) malloc(10*sizeof(char)); 
    if (f == 0)
      sprintf(geom_label[f],"GEOMETRY");
    else
      sprintf(geom_label[f],"GEOMETRY%d",f+1);

    num_elem = 0;
    ip_count(geom_label[f], &num_elem, 0);

    simple_geom = 1;
    entry_length = 0;
    for(i=0; i < num_elem; i++) {
      ip_count(geom_label[f], &entry_length, 1, i);
      if(entry_length > 1) simple_geom = 0;
    }

    if(simple_geom && num_elem%4) {
      sprintf(error_message, "Problem with number of elements in %s.",geom_label[f]);
      punt(error_message);
    }
    if(simple_geom) {
      frag_num_allatoms[f] = num_elem/4;
    }
    else {
      frag_num_allatoms[f] = 0;
      ip_count(geom_label[f],&(frag_num_allatoms[f]),0);
    }
    if (frag_num_allatoms[f] == 0) {
      sprintf(error_message, "%s is empty!",geom_label[f]);
      punt(error_message);
    }

    num_allatoms += frag_num_allatoms[f];

    /* Figure out how many non-dummy atoms are there */
    frag_num_atoms[f] = 0;
    for(i=0;i<frag_num_allatoms[f];i++) {
      if(simple_geom) 
        errcod = ip_string(geom_label[f], &atom_label, 1, i*4);
      else
        errcod = ip_string(geom_label[f], &atom_label, 2, i,0);
      if (errcod != IPE_OK) {
        sprintf(error_message, "Problem reading %s array!", geom_label[f]);
        punt(error_message);
      }
      if (strcmp(atom_label,"X"))
        ++frag_num_atoms[f];
    }
    num_atoms += frag_num_atoms[f];
  }

  if (num_allatoms > MAXATOM)
    punt("There are more atoms than allowed!");

  frag_atom[0] = frag_allatom[0] = 0;
  for (f=1; f<nfragments; ++f) {
    frag_atom[f] = frag_atom[f-1] + frag_num_atoms[f-1];
    frag_allatom[f] = frag_allatom[f-1] + frag_num_allatoms[f-1];
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
  all_atomcount = 0;
  for (f=0; f<nfragments; ++f) {

    for(i=0;i<frag_num_allatoms[f];i++){
      if(simple_geom) 
        errcod = ip_string(geom_label[f],&atom_label,1,4*i);
      else
        errcod = ip_string(geom_label[f],&atom_label,2,i,0);
      if (errcod != IPE_OK) {
        sprintf(error_message,"Problem with the %s array.", geom_label[f]);
        punt(error_message);
      }
      if (strcmp(atom_label,"X")) {
         atom_num(atom_label, &Z);
         free(atom_label);
         elemsymb_charges[atomcount] = Z;
         element[atomcount] = elem_name[(int)Z];
         full_element[all_atomcount] = elem_name[(int)Z];
         geometry[atomcount] = full_geom[all_atomcount];
         atom_dummy[all_atomcount] = 0;
         ++atomcount;
       }
       else {
         static const char* dummy_elem_name = "X";
         full_element[all_atomcount] = const_cast<char*>(dummy_elem_name);
         free(atom_label);
         atom_dummy[all_atomcount] = 1;
       }
  
      for(j=0; j<3;j++){
        if(simple_geom) 
          errcod = ip_data(geom_label[f],"%lf", &tmp,1,4*i+j+1);
        else
          errcod = ip_data(geom_label[f],"%lf", &tmp,2,i,j+1);
        if (errcod != IPE_OK) {
          sprintf(error_message,"Problem with the %s array.", geom_label[f]);
          punt(error_message);
        }
        else
          full_geom[all_atomcount][j] = tmp*conv_factor;
      }
      ++all_atomcount;
    }
  }

  for (f=0; f<nfragments; ++f)
    free(geom_label[f]);
  free(geom_label);
  read_charges();

  return;
}

}} // namespace psi::input
