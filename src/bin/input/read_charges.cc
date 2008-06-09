/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
/*
** This functions gets an array of user specified charges (if it exists)
** Added to facilitate counterpoise corrections with ghost atoms
** July-2001 GST
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

void read_charges()
{
  int i, j, errcod, atomcount, f;
  char charge_lbl[20], error_message[80];

  nuclear_charges = init_array(num_atoms);
// element gets set by read_zmat, read_cart or read_geomdat
// so there is no point in doing it here (RAK, 2008)
//  element = (char **) malloc(sizeof(char *)*num_atoms);

  atomcount = 0;

  for (f=0; f<nfragments; ++f) {

    if (f == 0)
      sprintf(charge_lbl,"CHARGES");
    else
      sprintf(charge_lbl,"CHARGES%d",f+1);

    if( ip_exist(charge_lbl,0) ) {
      ip_count(charge_lbl, &i, 0) ;
      if(i != frag_num_atoms[f]) {
        sprintf(error_message,"Number of charges not equal to number of atoms (excluding dummy) in %s!",charge_lbl);
        punt(error_message);
      }
      errcod = ip_double_array(charge_lbl, nuclear_charges+atomcount, frag_num_atoms[f]) ;
      if (errcod != IPE_OK) {
        sprintf(error_message,"Problem reading the %s array!", charge_lbl);
        punt(error_message);
      }
    }
    /* IF USER DOES NOT SPECIFY CHARGES, POINT TO DEFAULT CHARGES */
    else {
      for(i=0;i<frag_num_atoms[f];i++) {
        nuclear_charges[i+atomcount] = elemsymb_charges[i+atomcount];
//        element[i+atomcount] = elem_name[(int)elemsymb_charges[i+atomcount]];
      }
    }

    atomcount += frag_num_atoms[f];
  }
}

}} // namespace psi::input
