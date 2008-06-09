/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.9  2004/05/03 04:32:40  crawdad
 * Major mods based on merge with stable psi-3-2-1 release.  Note that this
 * version has not been fully tested and some scf-optn test cases do not run
 * correctly beccause of changes in mid-March 2004 to optking.
 * -TDC
 *
/* Revision 1.8.8.1  2004/04/10 19:41:32  crawdad
/* Fixed the DIIS code for UHF cases.  The new version uses the Pulay scheme of
/* building the error vector in the AO basis as FDS-SDF, followed by xformation
/* into the orthogonal AO basis.   This code converges faster for test cases
/* like cc8, but fails for linearly dependent basis sets for unknown reasons.
/* -TDC
/*
/* Revision 1.8  2002/04/03 02:06:01  janssen
/* Finish changes to use new include paths for libraries.
/*
/* Revision 1.7  2002/03/25 02:51:57  janssen
/* libciomr.h -> libciomr/libciomr.h
/*
/* Revision 1.6  2000/10/13 19:51:20  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.5  2000/07/10 18:03:31  sbrown
/* Enabling cscf to send over just the occupied SCF eigenvector for DFT
/* calculations.  Only done for the RHF case.
/*
/* Revision 1.4  2000/07/06 21:06:05  sbrown
/* Fixed a seg fault inf form_vec.c
/*
/* Revision 1.3  2000/07/06 20:04:01  sbrown
/* Added capabilities to send the eigenvector to cints for DFT
/* calculations.
/*
/* Revision 1.2  2000/07/05 21:47:30  sbrown
/* Enabled the code to export the SCF eigenvector to CINTS when doing DFT.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.2  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: form_vec.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

namespace psi { namespace cscf {

void form_vec()
{
  int i,nn,num_mo;
  int j,k,l;
  double **ctrans;
  double **temp;
  double **sqhmat;
  double **sqhmat2;
  double tol=1.0e-20;
  struct symm *s;
   

  ctrans = (double **) init_matrix(nsfmax,nsfmax); 
  temp = (double **) init_matrix(nsfmax,nsfmax);
  sqhmat = (double **) init_matrix(nsfmax,nsfmax);
  sqhmat2 = (double **) init_matrix(nsfmax,nsfmax);
  double* hevals = init_array(nsfmax);
   
  for (i=0; i < num_ir ; i++) {
    s = &scf_info[i];
    if (nn=s->num_so) {
      num_mo = s->num_mo;
      tri_to_sq(s->hmat,sqhmat,nn);
      mmult(s->sahalf,1,sqhmat,0,temp,0,num_mo,nn,nn,0);
      mmult(temp,0,s->sahalf,0,sqhmat2,0,num_mo,nn,num_mo,0);
      // This hack is borrowed from MPQC. Important enough to implement!
      if (hcore_guess == 1 && num_mo == nn) {
        // MPQC uses eigenvalues of the core hamiltonian in a nonorthogonal basis
        // to assign occupations. It seems to work better than the usual method.
        // right now this is an undocumented capability and not enabled by default
        sq_rsp(num_mo,num_mo,sqhmat,s->hevals,1,ctrans,tol);
        sq_rsp(num_mo,num_mo,sqhmat2,hevals,1,ctrans,tol);
      }
      else {
        sq_rsp(num_mo,num_mo,sqhmat2,s->hevals,1,ctrans,tol);
      }
      mxmb(s->sahalf,1,0,ctrans,1,0,s->cmat,1,0,nn,num_mo,num_mo);
      for(k=0;k < nn; k++)
        for(l=0;l < num_mo; l++) s->sahalf[k][l]=s->cmat[k][l];
      
      if (print & 2) {
        double* occs = init_array(num_mo);
        fprintf(outfile, "\nguess vector for irrep %s\n", s->irrep_label);
        //print_mat(s->cmat, nn, num_mo, outfile);
        eigout(s->cmat, s->hevals, occs, nn, num_mo, outfile);
        free(occs);
      }
    }
  }
  
  free(hevals);
  free_matrix(ctrans,nsfmax);
  free_matrix(temp,nsfmax);
  free_matrix(sqhmat,nsfmax);
  free_matrix(sqhmat2,nsfmax);
}

}} // namespace psi::cscf
