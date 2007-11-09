/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

/*
** PHASE.C
**
** This routine forms the product: C_new(t)*S*C_old = ~I to check the
** phases of the MOs.  When a phase change occurs, the new MO is
** corrected.  If MO swapping occurs, the phase_check flag is set to zero
** and later written to file30.  This prevents the correlated routines
** from trying to restart from the old wavefunction.
**
** T. Daniel Crawford 6/19/96
**
** Modified to include UHF orbitals.
** -TDC, 5/03
**
*/

int phase(void)
{
  int i, j, k, m;
  int nn, num_mo, row_max;
  int phase_chk;
  double maxvalue;
  double **smat, **tmp, **identity;
  double ***cnew;
  struct symm *s;

  phase_chk = 1;
  maxvalue = 0.0;

  cnew = (double ***) malloc(num_ir * sizeof(double **));

  if(uhf) {
    for(m=0; m < 2; m++) {
      for(k=0; k < num_ir; k++) {
	s = &scf_info[k];
	if(nn=s->num_so) {
	  num_mo = s->num_mo;
	  cnew[k] = block_matrix(nn, num_mo);
	  for(j=0; j < nn; j++)
	    for(i=0; i < num_mo; i++)
	      cnew[k][j][i] = spin_info[m].scf_spin[k].cmat[j][i];
	}
      }

      for(k=0; k < num_ir; k++) {
	s = &scf_info[k];
	if(nn=s->num_so) {
	  num_mo = s->num_mo;

	  smat = block_matrix(nn, nn);
	  tmp = block_matrix(num_mo, nn);
	  identity = block_matrix(num_mo, num_mo);

	  tri_to_sq(s->smat, smat, nn);

	  mmult(cnew[k],1,smat,0,tmp,0,num_mo,nn,nn,0);
	  mmult(tmp,0,spin_info[m].scf_spin[k].cmat_orig,0,identity,0,num_mo,nn,num_mo,0);

	  for(j=0; j < num_mo; j++) {
	    maxvalue = 0.0;
	    for(i=0; i < num_mo; i++) {
	      if(fabs(identity[j][i]) > maxvalue) {
		maxvalue = fabs(identity[j][i]);
		row_max = i;
	      }
	    }
	    if(row_max != j) phase_chk = 0;
	  }

	  if(phase_chk) {
	    for(i=0; i < num_mo; i++) {
	      if(identity[i][i] < 0.0)
		for(j=0; j < nn; j++) cnew[k][j][i] = -(cnew[k][j][i]);
	    }
	  }

	  free_block(smat);
	  free_block(tmp);
	  free_block(identity);

	}
      }

      if(phase_chk) {
	fprintf(outfile, "\n Correcting phases of orbitals of spin type %1d.\n", m);
	for(k=0; k < num_ir; k++) {
	  s = &scf_info[k];
	  if(nn=s->num_so) {
	    num_mo = s->num_mo;
	    for(j=0; j < nn; j++)
	      for(i=0; i < num_mo; i++)
		spin_info[m].scf_spin[k].cmat[j][i] = cnew[k][j][i];
	    free_block(cnew[k]);
	  }
	}
      }
      else fprintf(outfile, "\n No phase correction for spin type %1d possible.\n", m);

    }
  }
  else {
    /* Make a copy of the new MO vector */
    for(k=0; k < num_ir; k++) {
      s = &scf_info[k];
      if(nn=s->num_so) {
	num_mo = s->num_mo;
	cnew[k] = init_matrix(nn,num_mo);
	for(j=0; j < nn; j++)
	  for(i=0; i < num_mo; i++)
	    cnew[k][j][i] = s->cmat[j][i];
      }
      /*      fprintf(outfile, "MOs for Irrep %d\n", k);
	      print_mat(cnew[k], nn, num_mo, outfile); */
    }

    for(k=0; k < num_ir; k++) {
      s = &scf_info[k];
      if(nn=s->num_so) {
	num_mo = s->num_mo;

	smat = init_matrix(nn,nn);
	tmp = init_matrix(num_mo,nn);
	identity = init_matrix(num_mo,num_mo);

	/* Unpack the overlap matrix */
	tri_to_sq(s->smat,smat,nn);

	/* Form ~I = C^t(new) * S * C(old) */
	/*	  mxmb(cnew[k],nn,1,smat,1,nn,tmp,1,nn,nn,nn,nn);
		  mxmb(tmp,1,nn,s->cmat_orig,1,nn,identity,1,nn,nn,nn,nn);*/
	mmult(cnew[k],1,smat,0,tmp,0,num_mo,nn,nn,0);
	mmult(tmp,0,s->cmat_orig,0,identity,0,num_mo,nn,num_mo,0);

	/*	  fprintf(outfile, "Approximate Identity Matrix for Irrep %d\n", k);
		  print_mat(identity, num_mo, num_mo, outfile); */

	/* Check for MO swapping */
	for(j=0; j < num_mo; j++) {
	  maxvalue = 0.0;
	  for(i=0; i < num_mo; i++) {
	    if(fabs(identity[j][i]) > maxvalue) {
	      maxvalue = fabs(identity[j][i]);
	      row_max = i;
	    }
	  }
	  if(row_max != j) phase_chk = 0;
	}

	/* Now correct the MO phases, if necessary */
	if(phase_chk) {
	  for(i=0; i < num_mo; i++) {
	    if(identity[i][i] < 0.0) {
	      for(j=0; j < nn; j++) cnew[k][j][i] = -(cnew[k][j][i]);
	    }
	  }
	  /*	      fprintf(outfile, "Corrected MOs for irrep %d\n", k);
		      print_mat(cnew[k], nn, num_mo, outfile); */
	}
	free_matrix(smat,nn);
	free_matrix(tmp,num_mo);
	free_matrix(identity,num_mo);
      }
    }

    /* Finally, put the corrected MOs back into s->cmat if no swapping
       occurred */
    if(phase_chk) {
      fprintf(outfile, "\n Correcting phases of orbitals.\n");
      for(k=0; k < num_ir; k++) {
	s = &scf_info[k];
	if(nn=s->num_so) {
	  num_mo = s->num_mo;
	  for(j=0; j < nn; j++)
	    for(i=0; i < num_mo; i++)
	      s->cmat[j][i] = cnew[k][j][i];
	  free_matrix(cnew[k],nn);
	}
      }
    }
    else fprintf(outfile, "\n No phase correction possible.\n");

  }

  free(cnew);
  return(phase_chk);

}

}} // namespace psi::cscf
