/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.7  2004/10/08 17:33:41  nruss
 * Modified iwl_rdone() to acept boolean for deleteing one-electron ints. -NJR
 *
/* Revision 1.6  2002/12/06 15:50:32  crawdad
/* Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
/* necessary.  This is new for the PSI3 execution driver.
/* -TDC
/*
/* Revision 1.5  2002/04/03 02:06:01  janssen
/* Finish changes to use new include paths for libraries.
/*
/* Revision 1.4  2002/03/06 22:44:41  sherrill
/* Add new keyword orthog_only = true to just orthogonalize orbitals and do
/* nothing else.
/*
/* Revision 1.3  2001/06/21 21:00:37  crawdad
/* I have simplified the libiwl functions iwl_rdone() and iwl_wrtone() to only
/* read and write one-electron quantities and to more explicitly use the libpsio
/* structure to allow multiple quantities in a single one-electron IWL file.
/* The frozen-core energy is no longer dealt with in these functions, but is
/* now handled in libfile30.  The argument lists for these functions have
/* therefore changed quite a lot, and I've tried to correct all the PSI3
/* codes that are affected.
/* -TDC
/*
/* Revision 1.2  2000/10/13 19:51:21  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:31  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/11/02 18:10:14  evaleev
/* Direct SCF improved
/*
/* Revision 1.2  1999/08/17 19:04:16  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: rdone.c 3444 2007-08-03 20:49:20Z arnstein $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libiwl/iwl.h>

namespace psi { namespace cscf {

void rdone_iwl()
{
  int stat;
  int ntri = ioff[nbasis];
  int i,j,k,jj,kk;
  int max,off;
  double *ints;
  double e_fzc;

  /* If it's a direct SCF run - tell CINTS to compute one-electron integrals */
  if (direct_scf) {
    stat = system("cints --oeints");
    switch (stat) {
    case 0:
      /* CINTS ran successfully - continue */
      break;

    default:
      /* Something went wrong */
      fprintf(outfile,"  rdone_iwl: System call to CINTS failed. Check to see if it's in your PATH\n");
      fprintf(stderr,"System call to CINTS failed. Check to see if it's in your PATH.\n");
      exit(PSI_RETURN_FAILURE);
    }
  }

  ints = init_array(ntri);

  /* S integrals */
  stat = iwl_rdone(itapS,PSIF_SO_S,ints,ntri,0, 0, outfile); 
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].smat[ioff[j]+k] = ints[ioff[jj]+kk];
      }
    }
  }

  /* T integrals */
  stat = iwl_rdone(itapT,PSIF_SO_T,ints,ntri, 0, 0, outfile);
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].tmat[ioff[j]+k] = ints[ioff[jj]+kk];
      }
    }
  }

  /* V integrals */
  stat = iwl_rdone(itapV,PSIF_SO_V,ints,ntri,0, 0, outfile);
  for(i=0;i<num_ir;i++) {
    max = scf_info[i].num_so;
    off = scf_info[i].ideg;
    for(j=0;j<max;j++) {
      jj = j + off;
      for(k=0;k<=j;k++) {
	kk = k + off;
	scf_info[i].hmat[ioff[j]+k] = ints[ioff[jj]+kk] + scf_info[i].tmat[ioff[j]+k];
      }
    }
  }
  free(ints);
  ints = NULL;
}

}} // namespace psi::cscf
