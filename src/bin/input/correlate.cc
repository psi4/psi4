/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
/** This function returns a lookup array that correlates irreps of a point group
 ** with irreps of a subgroup produced by a displacement of a particular
 ** irrep (in the higher point group). It returns an array that contains the
 ** the irrep of the sugroup that corresponds with the irreps of the higher
 ** point group.
 ** returns: int *array with dimension (number of old irreps) */ 
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <libciomr/libciomr.h>
#include "global.h"

namespace psi { namespace input {

int *correlate(char *ptgrp, int irrep, int *nirreps_old, int *nirreps_new)
{
  int i;
  int *arr;

  if (strcmp(ptgrp,"C1 ") == 0)
    *nirreps_old = 1;
  else if (strcmp(ptgrp,"Cs ") == 0)
    *nirreps_old = 2;
  else if (strcmp(ptgrp,"Ci ") == 0)
    *nirreps_old = 2;
  else if (strcmp(ptgrp,"C2 ") == 0)
    *nirreps_old = 2;
  else if (strcmp(ptgrp,"C2v") == 0)
    *nirreps_old = 4;
  else if (strcmp(ptgrp,"D2 ") == 0)
    *nirreps_old = 4;
  else if (strcmp(ptgrp,"C2h") == 0)
    *nirreps_old = 4;
  else if (strcmp(ptgrp,"D2h") == 0)
    *nirreps_old = 8;
  else {
    fprintf(outfile,"point group %s unknown.\n",ptgrp);
    exit(1);
  }

  arr = init_int_array(*nirreps_old);

  if (irrep == 0) { /* return identity */
    *nirreps_new = *nirreps_old;
    for (i=0; i<*nirreps_old; ++i)
      arr[i] = i; 
    return arr;
  }

  *nirreps_new = *nirreps_old / 2;
  if ((strcmp(ptgrp,"C1 ") == 0) || (strcmp(ptgrp,"Cs ") == 0) ||
      (strcmp(ptgrp,"Ci ") == 0) || (strcmp(ptgrp,"C2 ") == 0) ) {
        arr[0] = 0;
  }
  else if ( (strcmp(ptgrp,"C2v") == 0) || (strcmp(ptgrp,"D2 ") == 0) ||
            (strcmp(ptgrp,"C2h") == 0) ) {
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
    }
  }
  else if (strcmp(ptgrp,"D2h") == 0) {
    /* 1,2,3 give C2h displaced geometries */
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 2; arr[6] = 3;  arr[7] = 3;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 3; arr[6] = 2;  arr[7] = 3;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
      arr[4] = 2;  arr[5] = 3; arr[6] = 3;  arr[7] = 2;
    }
    /* 4 gives D2 displaced geometries */
    else if (irrep == 4) { /* D2 */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 0;  arr[5] = 1; arr[6] = 2;  arr[7] = 3;
    }
    /* displacements along irreps 5,6,7 make C2v structures */
    /* care is taken to make sure definition of b1 and b2 will
       match those that input will generate - the following seems to work:
       b1u disp: has C2(z), b2 irrep symmetric wrt sigma(yz)
       b2u disp: has C2(y), b2 irrep symmetric wrt sigma(xy)
       b3u disp: has C2(x), b2 irrep symmetric wrt sigma(xz) */
    else if (irrep == 5) { /* b1u */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 1;  arr[5] = 0; arr[6] = 3;  arr[7] = 2;
    }
    else if (irrep == 6) { /* b2u */
      arr[0] = 0;  arr[1] = 3; arr[2] = 1;  arr[3] = 2;
      arr[4] = 1;  arr[5] = 2; arr[6] = 0;  arr[7] = 3;
    }
    else if (irrep == 7) { /* b3u */
      arr[0] = 0;  arr[1] = 2; arr[2] = 3;  arr[3] = 1;
      arr[4] = 1;  arr[5] = 3; arr[6] = 2;  arr[7] = 0;
    }
  }
  else {
    fprintf(outfile,"Point group unknown for correlation table.\n");
  }

  return arr;
}

}} // namespace psi::input
