/*!
  \file wrttwo.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include "iwl.h"

extern "C" {
	
/*!
** iwl_wrttwo()
**
** Write two electron ints to output in lexical order
** The "iwl" stands for "integrals with labels," and this is the proposed
** new standard for storing two-electron integrals and their (absolute)
** orbital labels.  This function closes the output file when finished.
**
**    \param itap     = unit to write to
**    \param nbfso    = number of basis functions in symmetry orbitals
**    \param ints     = two electron integrals 
**    \param ioff     = the old ioff array for lexical ordering
**    \param printflg = print flag (1 or 0)
**    \param outfile  =  output file
**
** Revised 6/27/96 by CDS
** \ingroup (IWL)
*/
void iwl_wrttwo(int itap, int nbfso, double *ints, int *ioff, double toler, 
                int printflg, FILE *outfile)
{
  struct iwlbuf Buf;

  iwl_buf_init(&Buf, itap, toler, 0, 0);
  iwl_buf_wrt_all(&Buf, nbfso, ints, ioff, printflg, outfile);
  iwl_buf_flush(&Buf, 1);
  iwl_buf_close(&Buf, 1);

}

} /* extern "C" */