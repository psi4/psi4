/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
/*
** This functions gets an array of user specified charges (if it exists)
** Added to facilitate counterpoise corrections with ghost atoms
** July-2001 GST
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

void read_charges()
{
  int i, j, errcod;

  /* INIT GLOBAL ARRAYS */
  nuclear_charges = init_array(num_atoms);
  element = (char **) malloc(sizeof(char *)*num_atoms);

  if( ip_exist("CHARGES",0) ) {
    ip_count("CHARGES", &i, 0) ;
    if(i != num_atoms) {
      punt("Number of charges not equal to number of atoms (excluding dummy)");
      }
    errcod = ip_double_array("CHARGES", nuclear_charges, num_atoms) ;
    if (errcod != IPE_OK) {
      punt("Problem reading the CHARGES array.");
      }
    for(i=0;i<num_atoms;i++)
      element[i] = elem_name[(int)nuclear_charges[i]];
  }
  /* IF USER DOES NOT SPECIFY CHARGES, POINT TO DEFAULT CHARGES */
  else {
    for(i=0;i<num_atoms;i++) {
      nuclear_charges[i] = elemsymb_charges[i];
      element[i] = elem_name[(int)elemsymb_charges[i]];
    }
  }
}

}} // namespace psi::input
