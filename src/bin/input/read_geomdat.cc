/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void read_geomdat()
{
  FILE *geomdat;
  int i, j, errcod;
  double Z = 0.0;
  double tmp = 0.0;
  char entry_name[20];

  ffile(&geomdat, "geom.dat", 2);
  if (geomdat != NULL) {
    ip_append(geomdat, outfile);
    fclose(geomdat);
  }

  num_allatoms = 0;
  sprintf(entry_name,"GEOMETRY%d",geomdat_entry);
  ip_count(entry_name,&num_allatoms,0);
  if (num_allatoms == 0)
    punt("The entry in geom.dat is empty or missing!");
  else if (num_allatoms > MAXATOM)
    punt("There are more atoms than allowed!");
  num_atoms = num_allatoms;

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  full_geom = block_matrix(num_allatoms,3);
  geometry = (double **) malloc(num_atoms*sizeof(double *));
  atom_dummy = (int *) malloc(sizeof(int)*num_allatoms);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  full_element = (char **) malloc(sizeof(char *)*num_allatoms);
  elemsymb_charges = init_array(num_atoms);

  for(i=0;i<num_atoms;i++){
    errcod = ip_data(entry_name,"%lf",&Z,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem with the geom.dat entry.");
    elemsymb_charges[i] = Z;
    //element[i] = elem_name[(int)Z];
    //full_element[i] = elem_name[(int)Z];
    for(j=0; j<3;j++){
      errcod = ip_data(entry_name,"%lf", &tmp,2,i,j+1);
      if (errcod != IPE_OK)
	punt("Problem with the geom.dat entry.");
      else
	full_geom[i][j] = tmp;
    }
    geometry[i] = full_geom[i];
    atom_dummy[i] = 0;
  }

  // --read_geomdat implies only 1 frament - but we create these here so
  // that read_charges() works and a general free call can be made later
  frag_num_atoms = (int *) malloc(nfragments*sizeof(int));
  frag_num_allatoms = (int *) malloc(nfragments*sizeof(int));
  frag_atom = (int *) malloc(nfragments*sizeof(int));
  frag_allatom = (int *) malloc(nfragments*sizeof(int));
  frag_num_atoms[0] = num_atoms;
  frag_num_allatoms[0] = num_allatoms;
  frag_atom[0] = 0;
  frag_allatom[0] = 0;

  read_charges();

  return;
}

}} // namespace psi::input
