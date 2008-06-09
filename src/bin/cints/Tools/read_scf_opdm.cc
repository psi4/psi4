/*! \file read_scf_opdm.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libint/libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"

#define OFFDIAG_DENS_FACTOR 0.5

namespace psi { namespace CINTS {

void read_scf_opdm()
{ 
  int i, j, ij, nao_i, nao_j, sh_i, sh_j;
  int ioffset, joffset;
  int nstri = ioff[Symmetry.num_so];
  double *dens, *denso;
  double **sq_dens, **sq_denso;
  double **tmp_mat;
  double pfac;
  double max_elem, temp;
  struct shell_pair* sp;

  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_read_entry(IOUnits.itapDSCF, "Integrals cutoff", (char *) &(UserOptions.cutoff), sizeof(double));
  dens = init_array(nstri);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    denso = init_array(nstri);

  if (UserOptions.reftype == uhf) {
      denso = init_array(nstri);
      psio_read_entry(IOUnits.itapDSCF, "Difference Alpha Density", (char *) dens, sizeof(double)*nstri);
      psio_read_entry(IOUnits.itapDSCF, "Difference Beta Density", (char *) denso, sizeof(double)*nstri);
      /*--- Form total and open-shell (spin) density matrices first ---*/
      for(i=0;i<nstri;i++) {
	  temp = dens[i] + denso[i];
	denso[i] = dens[i] - denso[i];
	dens[i] = temp;
      }
  }
  else {
      psio_read_entry(IOUnits.itapDSCF, "Difference Density", (char *) dens, sizeof(double)*nstri);
      if (UserOptions.reftype == rohf) {
	  denso = init_array(nstri);
	  psio_read_entry(IOUnits.itapDSCF, "Difference Open-Shell Density", (char *) denso, sizeof(double)*nstri);
      }
  }
  psio_close(IOUnits.itapDSCF, 1);

  
  /*--------------------------
    Remove after done testing
   --------------------------*/
/*  fprintf(outfile,"  Total density matrix in SO basis :\n");
  print_array(dens,Symmetry.num_so,outfile);
  fprintf(outfile,"\n\n");
  if (UserOptions.reftype == rohf) {
    fprintf(outfile,"  Total open-shell density matrix in SO basis :\n");
    print_array(denso,Symmetry.num_so,outfile);
    fprintf(outfile,"\n\n");
  }*/

  
  /*------------------------
    convert to square forms
   ------------------------*/
  sq_dens = init_matrix(BasisSet.num_ao,BasisSet.num_ao);
  ij = 0;
  for(i=0;i<Symmetry.num_so;i++) {
    for(j=0;j<i;j++) {
      sq_dens[i][j] = sq_dens[j][i] =  OFFDIAG_DENS_FACTOR*dens[ij];
      ij++;
    }
    sq_dens[i][i] = dens[ij];
    ij++;
  }
  free(dens);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    sq_denso = init_matrix(BasisSet.num_ao,BasisSet.num_ao);
    ij = 0;
    for(i=0;i<Symmetry.num_so;i++) {
      for(j=0;j<i;j++) {
	sq_denso[i][j] = sq_denso[j][i] = OFFDIAG_DENS_FACTOR*denso[ij];
	ij++;
      }
      sq_denso[i][i] = denso[ij];
      ij++;
    }
    free(denso);
  }

  /*----------------------
    transform to AO basis
   ----------------------*/
  tmp_mat = init_matrix(Symmetry.num_so,BasisSet.num_ao);
  mmult(sq_dens,0,Symmetry.usotao,0,tmp_mat,0,Symmetry.num_so,Symmetry.num_so,BasisSet.num_ao,0);
  mmult(Symmetry.usotao,1,tmp_mat,0,sq_dens,0,BasisSet.num_ao,Symmetry.num_so,BasisSet.num_ao,0);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    mmult(sq_denso,0,Symmetry.usotao,0,tmp_mat,0,Symmetry.num_so,Symmetry.num_so,BasisSet.num_ao,0);
    mmult(Symmetry.usotao,1,tmp_mat,0,sq_denso,0,BasisSet.num_ao,Symmetry.num_so,BasisSet.num_ao,0);
  }
  free_matrix(tmp_mat,Symmetry.num_so);


  /*--------------------------
    Remove after done testing
   --------------------------*/
  /*fprintf(outfile,"  Total density matrix in AO basis :\n");
  print_mat(sq_dens,BasisSet.num_ao,BasisSet.num_ao,outfile);
  fprintf(outfile,"\n\n");
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    fprintf(outfile,"  Total open-shell density matrix in AO basis :\n");
    print_mat(sq_denso,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n\n");
    }*/
  
  /*--------------------------------
    transform to shell-blocked form
   --------------------------------*/
  switch (UserOptions.reftype) {
    case rohf:
      for(sh_i=0;sh_i<BasisSet.num_shells;sh_i++) {
	ioffset = BasisSet.shells[sh_i].fao-1;
	nao_i = ioff[BasisSet.shells[sh_i].am];
	for(sh_j=0;sh_j<BasisSet.num_shells;sh_j++) {
	  sp = &(BasisSet.shell_pairs[sh_i][sh_j]);
	  joffset = BasisSet.shells[sh_j].fao-1;
	  nao_j = ioff[BasisSet.shells[sh_j].am];
	  if (sp->dmato == NULL)
	    sp->dmato = block_matrix(nao_i,nao_j);

	  max_elem = sp->Dmax;
	  for(i=0;i<nao_i;i++)
	    for(j=0;j<nao_j;j++) {
	      temp = sp->dmato[i][j] = sq_denso[ioffset+i][joffset+j];
	      temp = fabs(temp);
	      if (max_elem < temp)
		max_elem = temp;
	    }
	  sp->Dmax = max_elem;
	}
      }
      free_matrix(sq_denso,BasisSet.num_ao);

      for(sh_i=0;sh_i<BasisSet.num_shells;sh_i++) {
	ioffset = BasisSet.shells[sh_i].fao-1;
	nao_i = ioff[BasisSet.shells[sh_i].am];
	for(sh_j=0;sh_j<BasisSet.num_shells;sh_j++) {
	  sp = &(BasisSet.shell_pairs[sh_i][sh_j]);
	  joffset = BasisSet.shells[sh_j].fao-1;
	  nao_j = ioff[BasisSet.shells[sh_j].am];
	  if (sp->dmat == NULL)
	    sp->dmat = block_matrix(nao_i,nao_j);
	  max_elem = -1.0;
	  for(i=0;i<nao_i;i++)
	    for(j=0;j<nao_j;j++) {
	      temp = sp->dmat[i][j] = sq_dens[ioffset+i][joffset+j];
	      temp = fabs(temp);
	      if (max_elem < temp)
		max_elem = temp;
	    }
	  sp->Dmax = MAX(sp->Dmax,max_elem);
	}
      }
      free_matrix(sq_dens,BasisSet.num_ao);
      break;

    case uhf:
      for(sh_i=0;sh_i<BasisSet.num_shells;sh_i++) {
	ioffset = BasisSet.shells[sh_i].fao-1;
	nao_i = ioff[BasisSet.shells[sh_i].am];
	for(sh_j=0;sh_j<BasisSet.num_shells;sh_j++) {
	  sp = &(BasisSet.shell_pairs[sh_i][sh_j]);
	  joffset = BasisSet.shells[sh_j].fao-1;
	  nao_j = ioff[BasisSet.shells[sh_j].am];
	  /*--- Check if can reuse space from dmata ---*/
	  if (sp->dmato == NULL)
	      if (sp->dmata == NULL)
		  sp->dmato = block_matrix(nao_i,nao_j);
	      else {
		  sp->dmato = sp->dmata;
		  sp->dmata = NULL;
	      }
	  
	  max_elem = sp->Dmax;
	  for(i=0;i<nao_i;i++)
	      for(j=0;j<nao_j;j++) {
		  temp = sp->dmato[i][j] = sq_denso[ioffset+i][joffset+j];
		  temp = fabs(temp);
		  if (max_elem < temp)
		      max_elem = temp;
	      }
	  sp->Dmax = max_elem;
	}
      }
      free_matrix(sq_denso,BasisSet.num_ao);
      
      for(sh_i=0;sh_i<BasisSet.num_shells;sh_i++) {
	  ioffset = BasisSet.shells[sh_i].fao-1;
	  nao_i = ioff[BasisSet.shells[sh_i].am];
	  for(sh_j=0;sh_j<BasisSet.num_shells;sh_j++) {
	      sp = &(BasisSet.shell_pairs[sh_i][sh_j]);
	      joffset = BasisSet.shells[sh_j].fao-1;
	      nao_j = ioff[BasisSet.shells[sh_j].am];
	      /*--- Check if can reuse space from dmatb ---*/
	      if (sp->dmat == NULL)
		  if (sp->dmatb == NULL)
		      sp->dmat = block_matrix(nao_i,nao_j);
		  else {
		      sp->dmat = sp->dmatb;
		      sp->dmatb = NULL;
		  }
	      
	      max_elem = sp->Dmax;
	      for(i=0;i<nao_i;i++)
		  for(j=0;j<nao_j;j++) {
		      temp = sp->dmat[i][j] = sq_dens[ioffset+i][joffset+j];
		      temp = fabs(temp);
		      if (max_elem < temp)
			  max_elem = temp;
		  }
	      sp->Dmax = MAX(sp->Dmax,max_elem);      
	  }
      }
  
      free_matrix(sq_dens,BasisSet.num_ao);
      break;

    case rhf:
      for(sh_i=0;sh_i<BasisSet.num_shells;sh_i++) {
	ioffset = BasisSet.shells[sh_i].fao-1;
	nao_i = ioff[BasisSet.shells[sh_i].am];
	for(sh_j=0;sh_j<BasisSet.num_shells;sh_j++) {
	  sp = &(BasisSet.shell_pairs[sh_i][sh_j]);
	  joffset = BasisSet.shells[sh_j].fao-1;
	  nao_j = ioff[BasisSet.shells[sh_j].am];
	  if (sp->dmat == NULL)
	    sp->dmat = block_matrix(nao_i,nao_j);
	  max_elem = -1.0;
	  for(i=0;i<nao_i;i++)
	    for(j=0;j<nao_j;j++) {
	      temp = sp->dmat[i][j] = sq_dens[ioffset+i][joffset+j];
	      temp = fabs(temp);
	      if (max_elem < temp)
		max_elem = temp;
	    }
	  sp->Dmax = max_elem;
	}
      }
      free_matrix(sq_dens,BasisSet.num_ao);
      break;
  }

  return;
}

};};
