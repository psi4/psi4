/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void rhf_fock_build(double **fock, double **D);
void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

/* rotate(): Rotate the orbitals using a linear transformation matrix
** built from converged T1 amplitudes.  I still need to add spin-restricted
** Brueckner rotations from my 1997 paper.
**
** TDC, 5/03
*/

int rotate(void)
{
  int i, a, ii, aa, j, ij, b, p, q, I, J, A, B;
  int h, nirreps, nso, nmo, ntri, stat;
  dpdfile2 T1;
  double **U, **S, **X, *scratch;
  double *evals, *work, **SO_S, **MO_S;
  double **scf, **scf_new, **scf_a, **scf_b;
  double **scf_orig, **scf_a_orig, **scf_b_orig;
  double max;
  double **D, **D_a, **D_b; /* SCF densities */
  double **fock, **fock_a, **fock_b; /* Fock matrices (SO or MO basis) */
  double ***Foo, ***Fvv; /* occ-occ and vir-vir block of Fock matrix */
  int *offset;
  int phase_ok=1, max_col;

  chkpt_init(PSIO_OPEN_OLD);

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  offset = init_int_array(nirreps);
  for(h=1; h < nirreps; h++)
    offset[h] = offset[h-1] + moinfo.orbspi[h-1];

  /* First check to see if we've already converged the orbitals */
  max = 0.0;
  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.occpi[h]; i++)
	for(a=0; a < moinfo.virtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = fabs(T1.matrix[h][i][a]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.aoccpi[h]; i++)
	for(a=0; a < moinfo.avirtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = fabs(T1.matrix[h][i][a]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.boccpi[h]; i++)
	for(a=0; a < moinfo.bvirtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = fabs(T1.matrix[h][i][a]);

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
  }

  if(fabs(max) <= params.bconv) {
    fprintf(outfile, "\tBrueckner orbitals converged.  Maximum T1 = %15.12f\n",
	    fabs(max));
    chkpt_close();
    return(1);
  }
  else 
    fprintf(outfile, "\tRotating orbitals.  Maximum T1 = %15.12f\n", fabs(max));

  /* grab the SO-basis overlap integrals for later use */
  SO_S = block_matrix(nso, nso);
  ntri = nso * (nso+1)/2;
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_S, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      SO_S[i][j] = SO_S[j][i] = scratch[ij];
    }
  free(scratch);

  if(params.ref == 0) { /* RHF */

    U = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) U[i][i] = 1.0;

    max = 0.0;
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.occpi[h]; i++) {
	ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[h]];
	for(a=0; a < moinfo.virtpi[h]; a++) {
	  aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    scf = chkpt_rd_scf();
    scf_orig = chkpt_rd_scf();
    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);
    free_block(U);
    free_block(scf);

    /* transform the overlap into the new MO basis */
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    /* build S^-1/2 for this basis */
    evals = init_array(nmo);
    work = init_array(nmo*3);
    if((stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3))) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }
    free(evals);
    free(work);
    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the new MO basis */
    scf = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf[0][0]),nmo);
    free_block(S);
    free_block(MO_S);
    free_block(scf_new);

    /* build the SO-basis density for the new MOs */
    D = block_matrix(nso,nso);
    for(h=0; h < nirreps; h++)
      for(p=offset[h]; p < offset[h]+moinfo.orbspi[h]; p++)
	for(q=offset[h]; q < offset[h]+moinfo.orbspi[h]; q++)
	  for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]+moinfo.occpi[h]; i++)
	    D[p][q] += scf[p][i] * scf[q][i];

    /* build the SO-basis Fock matrix */
    fock = block_matrix(nso, nso);
    rhf_fock_build(fock, D);
    free_block(D);

    /*
    fprintf(outfile, "\n\tSO-basis Fock matrix:\n");
    mat_print(fock, nso, nso, outfile);
    */

    /* transform the fock matrix to the new MO basis */
    X = block_matrix(nso,nso);
    C_DGEMM('n','n',nso,nmo,nso,1.0,&(fock[0][0]),nso,&(scf[0][0]),nmo,
	    0,&(X[0][0]),nso);
    C_DGEMM('t','n',nmo,nmo,nso,1.0,&(scf[0][0]),nmo,&(X[0][0]),nso,
	    0,&(fock[0][0]),nso);
    free_block(X);

    /*
    fprintf(outfile, "\n\tMO-basis Fock matrix:\n");
    mat_print(fock, nmo, nmo, outfile);
    */

    /* extract the occ-occ and vir-vir block of the Fock matrix */
    Foo = (double ***) malloc(nirreps * sizeof(double **));
    Fvv = (double ***) malloc(nirreps * sizeof(double **));
    X = block_matrix(nmo, nmo);
    for(h=0; h < nirreps; h++) {

      /* leave the frozen core orbitals alone */
      for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++) X[i][i] = 1.0;

      Foo[h] = block_matrix(moinfo.occpi[h], moinfo.occpi[h]);
      Fvv[h] = block_matrix(moinfo.virtpi[h], moinfo.virtpi[h]);

      for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.occpi[h]; i++,I++)
	for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.occpi[h]; j++,J++)
	  Foo[h][I][J] = fock[i][j];

      for(a=offset[h]+moinfo.frdocc[h]+moinfo.occpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	for(b=offset[h]+moinfo.frdocc[h]+moinfo.occpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	  Fvv[h][A][B] = fock[a][b];

      /*
      fprintf(outfile, "\n\tOcc-occ Fock matrix for irrep %d:\n", h);
      mat_print(Foo[h], moinfo.occpi[h], moinfo.occpi[h], outfile);

      fprintf(outfile, "\n\tVir-vir Fock matrix for irrep %d:\n", h);
      mat_print(Fvv[h], moinfo.virtpi[h], moinfo.virtpi[h], outfile);
      */

      if(moinfo.occpi[h]) {
	evals = init_array(moinfo.occpi[h]);
	work = init_array(3*moinfo.occpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.occpi[h], &(Foo[h][0][0]), 
			  moinfo.occpi[h], evals, work, moinfo.occpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in Foo[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	/*
	fprintf(outfile, "\n\tEigenfunctions of Occ-occ Fock matrix for irrep %1d:\n", h);
	mat_print(Foo[h], moinfo.occpi[h], moinfo.occpi[h], outfile);
	*/

	for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.occpi[h]; i++,I++)
	  for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.occpi[h]; j++,J++)
	    X[i][j] = Foo[h][J][I];
      }

      if(moinfo.virtpi[h]) {
	evals = init_array(moinfo.virtpi[h]);
	work = init_array(3*moinfo.virtpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.virtpi[h], &(Fvv[h][0][0]), moinfo.virtpi[h], 
			  evals, work, moinfo.virtpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in Fvv[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	/*
	fprintf(outfile, "\n\tEigenfunctions of Vir-vir Fock matrix for irrep %1d:\n", h);
	mat_print(Fvv[h], moinfo.virtpi[h], moinfo.virtpi[h], outfile);
	*/

	for(a=offset[h]+moinfo.frdocc[h]+moinfo.occpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	  for(b=offset[h]+moinfo.frdocc[h]+moinfo.occpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	    X[a][b] = Fvv[h][B][A];
      }

      free_block(Foo[h]);
      free_block(Fvv[h]);
    }
    free(Foo);
    free(Fvv);
    free_block(fock);

    /* semicanonicalization of the basis */
    /*
    fprintf(outfile, "\n\tSemicanonical transformation matrix:\n");
    mat_print(X, nmo, nmo, outfile);
    */

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','n', nso, nmo, nmo, 1, &(scf[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(scf_new[0][0]), nmo);
    free_block(X);
    free_block(scf);

    /* Reorder new MO's to Pitzer and write to chkpt */
    /*
    fprintf(outfile, "\n\tSemicanonical Brueckner orbitals (Pitzer order):\n");
    mat_print(scf_new, nso, nmo, outfile);
    */

    /* correct orbital phases for amplitude restarts */
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_orig[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    scf_new[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);

    /*
    fprintf(outfile, "\n\tOriginal SCF MOs:\n");
    mat_print(scf_orig, nso, nmo, outfile);
    fprintf(outfile, "\n\tNew SCF MOs:\n");
    mat_print(scf_new, nso, nmo, outfile);
    */

    chkpt_wt_scf(scf_new);
    free_block(scf_new);
    free_block(scf_orig);

  }
  else if(params.ref == 2) { /* UHF */

    /* AA block */
    U = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) U[i][i] = 1.0;

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.aoccpi[h]; i++) {
	ii = moinfo.qt2pitzer_a[moinfo.qt_aocc[i] + moinfo.aocc_off[h]];
	for(a=0; a < moinfo.avirtpi[h]; a++) {
	  aa = moinfo.qt2pitzer_a[moinfo.qt_avir[a] + moinfo.avir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    scf = chkpt_rd_alpha_scf();
    scf_a_orig = chkpt_rd_alpha_scf();

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);
    free_block(U);
    free_block(scf);

    MO_S = block_matrix(nmo, nmo);

    /* transform the overlap into the new alpha MO basis */
    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    evals = init_array(nmo);
    work = init_array(nmo*3);
    if((stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3))) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }

    /* build S^-1/2 for this basis */
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }
    free(evals);
    free(work);
    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the basis */
    scf_a = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf_a[0][0]),nmo);
    free_block(S);
    free_block(MO_S);
    free_block(scf_new);

    /* BB block */
    U = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) U[i][i] = 1.0;

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.boccpi[h]; i++) {
	ii = moinfo.qt2pitzer_b[moinfo.qt_bocc[i] + moinfo.bocc_off[h]];
	for(a=0; a < moinfo.bvirtpi[h]; a++) {
	  aa = moinfo.qt2pitzer_b[moinfo.qt_bvir[a] + moinfo.bvir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    scf = chkpt_rd_beta_scf();
    scf_b_orig = chkpt_rd_beta_scf();

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);
    free_block(U);
    free_block(scf);

    MO_S = block_matrix(nmo, nmo);

    /* transform the overlap into the new beta MO basis */
    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    evals = init_array(nmo);
    work = init_array(nmo*3);
    if((stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3))) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }

    /* build S^-1/2 for this basis */
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }
    free(evals);
    free(work);
    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the basis */
    scf_b = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf_b[0][0]),nmo);
    free_block(S);
    free_block(MO_S);
    free_block(scf_new);

    /* build the SO-basis alpha and beta densities for the new MOs */
    D_a = block_matrix(nso, nso);
    D_b = block_matrix(nso, nso);
    for(h=0; h < nirreps; h++)
      for(p=offset[h]; p < offset[h]+moinfo.orbspi[h]; p++)
	for(q=offset[h]; q < offset[h]+moinfo.orbspi[h]; q++) {
	  for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h]; i++)
	    D_a[p][q] += scf_a[p][i] * scf_a[q][i];
	  for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h]; i++)
	    D_b[p][q] += scf_b[p][i] * scf_b[q][i];
	}

    /* build the alpha and beta SO-basis Fock matrices */
    fock_a = block_matrix(nso, nso);
    fock_b = block_matrix(nso, nso);
    uhf_fock_build(fock_a, fock_b, D_a, D_b);
    free_block(D_a);
    free_block(D_b);

    /* transform the fock matrices to the new alpha and beta MO bases */
    X = block_matrix(nso,nso);
    C_DGEMM('n','n',nso,nmo,nso,1.0,&(fock_a[0][0]),nso,&(scf_a[0][0]),nmo,
	    0,&(X[0][0]),nso);
    C_DGEMM('t','n',nmo,nmo,nso,1.0,&(scf_a[0][0]),nmo,&(X[0][0]),nso,
	    0,&(fock_a[0][0]),nso);

    C_DGEMM('n','n',nso,nmo,nso,1.0,&(fock_b[0][0]),nso,&(scf_b[0][0]),nmo,
	    0,&(X[0][0]),nso);
    C_DGEMM('t','n',nmo,nmo,nso,1.0,&(scf_b[0][0]),nmo,&(X[0][0]),nso,
	    0,&(fock_b[0][0]),nso);
    free_block(X);

    /** alpha Fock semicanonicalization **/

    Foo = (double ***) malloc(nirreps * sizeof(double **));
    Fvv = (double ***) malloc(nirreps * sizeof(double **));
    X = block_matrix(nmo, nmo);
    for(h=0; h < nirreps; h++) {
      /* leave the frozen core orbitals alone */
      for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++) X[i][i] = 1.0;

      Foo[h] = block_matrix(moinfo.aoccpi[h], moinfo.aoccpi[h]);
      Fvv[h] = block_matrix(moinfo.avirtpi[h], moinfo.avirtpi[h]);

      for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h]; i++,I++)
	for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h]; j++,J++)
	  Foo[h][I][J] = fock_a[i][j];

      for(a=offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	for(b=offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	  Fvv[h][A][B] = fock_a[a][b];

      if(moinfo.aoccpi[h]) {
	evals = init_array(moinfo.aoccpi[h]);
	work = init_array(3*moinfo.aoccpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.aoccpi[h], &(Foo[h][0][0]), 
			  moinfo.aoccpi[h], evals, work, moinfo.aoccpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h]; i++,I++)
	  for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h]; j++,J++)
	    X[i][j] = Foo[h][J][I];
      }

      if(moinfo.avirtpi[h]) {
	evals = init_array(moinfo.avirtpi[h]);
	work = init_array(3*moinfo.avirtpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.avirtpi[h], &(Fvv[h][0][0]), moinfo.avirtpi[h], 
			  evals, work, moinfo.avirtpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	for(a=offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	  for(b=offset[h]+moinfo.frdocc[h]+moinfo.aoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	    X[a][b] = Fvv[h][B][A];
      }

      free_block(Foo[h]);
      free_block(Fvv[h]);
    }
    free(Foo);
    free(Fvv);
    free_block(fock_a);

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','n', nso, nmo, nmo, 1, &(scf_a[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(scf_new[0][0]), nmo);
    free_block(X);
    free_block(scf_a);

    /* correct orbital phases for amplitude restarts */
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_a_orig[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    scf_new[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);

    chkpt_wt_alpha_scf(scf_new);
    free_block(scf_new);
    free_block(scf_a_orig);

    /** beta Fock semicanonicalization **/

    Foo = (double ***) malloc(nirreps * sizeof(double **));
    Fvv = (double ***) malloc(nirreps * sizeof(double **));
    X = block_matrix(nmo, nmo);
    for(h=0; h < nirreps; h++) {
      /* leave the frozen core orbitals alone */
      for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++) X[i][i] = 1.0;

      Foo[h] = block_matrix(moinfo.boccpi[h], moinfo.boccpi[h]);
      Fvv[h] = block_matrix(moinfo.bvirtpi[h], moinfo.bvirtpi[h]);

      for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h]; i++,I++)
	for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h]; j++,J++)
	  Foo[h][I][J] = fock_b[i][j];

      for(a=offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	for(b=offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	  Fvv[h][A][B] = fock_b[a][b];

      if(moinfo.boccpi[h]) {
	evals = init_array(moinfo.boccpi[h]);
	work = init_array(3*moinfo.boccpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.boccpi[h], &(Foo[h][0][0]), 
			  moinfo.boccpi[h], evals, work, moinfo.boccpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	for(i=offset[h]+moinfo.frdocc[h],I=0; i < offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h]; i++,I++)
	  for(j=offset[h]+moinfo.frdocc[h],J=0; j < offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h]; j++,J++)
	    X[i][j] = Foo[h][J][I];
      }

      if(moinfo.bvirtpi[h]) {
	evals = init_array(moinfo.bvirtpi[h]);
	work = init_array(3*moinfo.bvirtpi[h]);
	if((stat = C_DSYEV('v','u', moinfo.bvirtpi[h], &(Fvv[h][0][0]), moinfo.bvirtpi[h], 
			  evals, work, moinfo.bvirtpi[h]*3))) {
	  fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", 
		  h, stat);
	  exit(PSI_RETURN_FAILURE);
	}
	free(evals);
	free(work);

	for(a=offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
	  for(b=offset[h]+moinfo.frdocc[h]+moinfo.boccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
	    X[a][b] = Fvv[h][B][A];
      }

      free_block(Foo[h]);
      free_block(Fvv[h]);
    }
    free(Foo);
    free(Fvv);
    free_block(fock_b);

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','n', nso, nmo, nmo, 1, &(scf_b[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(scf_new[0][0]), nmo);
    free_block(X);
    free_block(scf_b);

    /* correct orbital phases for amplitude restarts */
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_b_orig[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    scf_new[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);

    chkpt_wt_beta_scf(scf_new);
    free_block(scf_new);
    free_block(scf_b_orig);

  }

  free_block(SO_S);
  free(offset);

  chkpt_wt_phase_check(phase_ok);
  chkpt_close();

  return 0;
}
}} // namespace psi::ccenergy
