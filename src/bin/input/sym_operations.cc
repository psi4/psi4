/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <cstdlib>
#include <cmath>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - C2_Z              */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/
void c2z(double *A, double *B)
{ 
  B[0] = -A[0];
  B[1] = -A[1];
  B[2] =  A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - C2_Y              */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/



void c2y(double *A, double *B)
{ 
  B[0] = -A[0];
  B[1] =  A[1];
  B[2] = -A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - C2_X              */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/



void c2x(double *A, double *B)
{ 
  B[0] =  A[0];
  B[1] = -A[1];
  B[2] = -A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - inversion         */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/



void inversion(double *A, double *B)
{ 
  B[0] = -A[0];
  B[1] = -A[1];
  B[2] = -A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - Sig_XY            */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/



void sig_xy(double *A, double *B)
{ 
  B[0] =  A[0];
  B[1] =  A[1];
  B[2] = -A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - Sig_YZ            */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/



void sig_yz(double *A, double *B)
{ 
  B[0] = -A[0];
  B[1] =  A[1];
  B[2] =  A[2];

  return;
}

/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - Sig_XZ            */
/*   which will be used to generate the non-symmetry unique atoms from the    */
/*   cartesians input by the user                                             */
/******************************************************************************/


void sig_xz(double *A, double *B)
{
  B[0] =  A[0];
  B[1] = -A[1];
  B[2] =  A[2];

  return;
}


/******************************************************************************/
/*   This function contains the Symmetry Operation Matrix - C4_Y              */
/*   which will be used to swap X and Z axis in reorient()                    */
/******************************************************************************/
void c4y(double *A, double *B)
{ 

   /* Perform C4 Operation on the old_coord. vector */
   B[0] = -A[2];
   B[1] =  A[1];
   B[2] =  A[0];

}

}} // namespace psi::input
