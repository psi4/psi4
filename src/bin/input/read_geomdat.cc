/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
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
    element[i] = elem_name[(int)Z];
    full_element[i] = elem_name[(int)Z];
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

  read_charges();

  return;
}

}} // namespace psi::input
