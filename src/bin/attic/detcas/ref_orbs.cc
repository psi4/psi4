/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** REF_ORBS.C
** 
** This file contains routines pertaining to the set of reference orbitals,
** C_0, from which the orbital rotation angles are defined.
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

#include <cstdlib>
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"


namespace psi { namespace detcas {

/*
** read_ref_orbs()
**
** This function reads in the ``reference'' orbitals, C_0, from a special
**  file.  If this file does not exist, then it is created, and the
**  current orbitals in file30 are placed there, and the array of orbital
**  rotation angles is reset to zero.
**
** Returns: 1 if read is successful, otherwise 0
*/
int read_ref_orbs(void)
{
  FILE *fp;
  int h, ir_orbs;

  ffileb_noexit(&fp,"orbs.dat",2);
  if (fp == NULL) {
    if (Params.print_lvl) 
      fprintf(outfile, "No orbs.dat file ... using new reference orbitals\n");
    return(0);
  }

  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_orbs = CalcInfo.orbs_per_irr[h];
    if (ir_orbs == 0) continue;
    if (fread(CalcInfo.mo_coeffs[h][0], sizeof(double), ir_orbs * ir_orbs,
              fp) != ir_orbs * ir_orbs) {
      fprintf(outfile, "Error reading reference orbitals.\n");
      fclose(fp);
      return(0);
    }
  }
    
  fclose(fp);
  return(1);
}


/*
** write_ref_orbs()
**
** This function initializes the set of reference orbitals.  The current
**  orbitals in file 30 are read and written to the reference orbital
**  disk file.
**
** Returns: 1 if successful, otherwise 0
*/
int write_ref_orbs(void)
{
  FILE *fp;
  int h, ir_orbs;

  ffileb_noexit(&fp,"orbs.dat",0);
  if (fp == NULL) {
    if (Params.print_lvl) 
      fprintf(outfile, "Can't open orbs.dat file!\n");
    return(0);
  }

  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_orbs = CalcInfo.orbs_per_irr[h];
    if (ir_orbs == 0) continue;
    if (fwrite(CalcInfo.mo_coeffs[h][0], sizeof(double), ir_orbs * ir_orbs,
              fp) != ir_orbs * ir_orbs) {
      fprintf(outfile, "Error writing reference orbitals.\n");
      fclose(fp);
      return(0);
    }
  }
    
  fclose(fp);
  return(1);
}

}} // end namespace psi::detcas

