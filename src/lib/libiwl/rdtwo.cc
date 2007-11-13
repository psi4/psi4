/*!
  \file rdtwo.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include "iwl.hpp"
#include "iwl.h"

  using namespace psi;
  
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

void IWL::rdtwo(PSIO* psio, int itap, double *ints, int *ioff, int norbs, 
      int nfzc, int nfzv, int printflg, FILE *outfile)
{
  IWL two_ints(psio, itap, 0.0, 1, 1);
  if ((nfzc == 0) && (nfzv == 0))
    two_ints.rd_all(ints, ioff, ioff, 0, ioff, printflg, outfile);
  else
    two_ints.rd_all_act(ints, ioff, ioff, 0, ioff, nfzc, norbs-nfzv-1, printflg, outfile);
  two_ints.set_keep(1);
}

extern "C" {	
/*!
** iwl_rdtwo(): read two electron ints from the given file.
** The "iwl" stands for "integrals with labels," and this is the proposed
** new standard for storing two-electron integrals and their (absolute)
** orbital labels.
**
**    \param itap     = unit to read from
**    \param ints     = two electron integrals (already allocated)
**    \param ioff     = the old ioff array for lexical ordering
**    \param norbs    = number of orbitals
**    \param nfzc     = number of frozen core orbitals
**    \param nfzv     = number of frozen virtual orbitals
**    \param printflg = print integrals as they're read 
**    \param outfile  = output file pointer
**
** David Sherrill, 1995
** \ingroup (IWL)
*/
void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs, 
      int nfzc, int nfzv, int printflg, FILE *outfile)
{
  struct iwlbuf Buf;
  
  iwl_buf_init(&Buf, itap, 0.0, 1, 1);
  if ((nfzc == 0) && (nfzv == 0))
    iwl_buf_rd_all(&Buf, ints, ioff, ioff, 0, ioff, printflg, outfile);
  else
    iwl_buf_rd_all_act(&Buf, ints, ioff, ioff, 0, ioff, nfzc, norbs-nfzv-1,
                       printflg, outfile);
  iwl_buf_close(&Buf, 1);
}

} /* extern "C" */