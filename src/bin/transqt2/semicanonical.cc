/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <exception.h>
#define EXTERN
#include "globals.h"

/*
** semicanonical_fock(): Compute the alpha- and beta-spin Fock
** matrices using the ROHF orbitals from the chkpt file and
** diagonalize the occ-occ and vir-vir blocks to form semicanonical
** orbitals, which are written in the checkpoint file for the
** subsequent transformation.  These orbitals are required for
** perturbation-based methods (e.g., (T), CC2, CC3, etc.) built upon
** an ROHF reference.
**
** This code is based on MLA's code in the old transqt, which was
** based on TDC's Brueckner code in ccenergy.
**
** TDC, 7/06
*/

namespace psi {
extern FILE* outfile;
  namespace transqt2 {

void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

void semicanonical_fock(void)
{
  int h, p, q, i, j, a, b, I, J, A, B;
  int stat, cnt;
  int *mo_offset, *so_offset, *aoccpi, *boccpi, *avirtpi, *bvirtpi;
  double **D_a, **D_b, **fock_a, **fock_b;
  double **X;
  double ***Foo, ***Fvv;
  double *evals, *work;
  double *alpha_evals, *beta_evals;
  double **C, **C_a, **C_b;
  int nirreps, nmo, nso;

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;

  aoccpi = init_int_array(nirreps);
  boccpi = init_int_array(nirreps);
  avirtpi = init_int_array(nirreps);
  bvirtpi = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
    aoccpi[h] = moinfo.clsdpi[h] + moinfo.openpi[h];
    boccpi[h] = moinfo.clsdpi[h];
    avirtpi[h] = moinfo.uoccpi[h];
    bvirtpi[h] = moinfo.uoccpi[h] + moinfo.openpi[h];
  }

  so_offset = init_int_array(nirreps);
  mo_offset = init_int_array(nirreps);
  for(h=1; h < nirreps; h++) mo_offset[h] = mo_offset[h-1] + moinfo.mopi[h-1];
  for(h=1; h < nirreps; h++) so_offset[h] = so_offset[h-1] + moinfo.sopi[h-1];

  /* build the SO-basis alpha and beta densities from the ROHF eigenvector */
  chkpt_init(PSIO_OPEN_OLD);
  C = chkpt_rd_scf();
  chkpt_close();

  D_a = block_matrix(nso, nso);
  D_b = block_matrix(nso, nso);
  for(h=0; h < nirreps; h++) {
    for(p=so_offset[h]; p < so_offset[h]+moinfo.sopi[h]; p++) {
      for(q=so_offset[h]; q < so_offset[h]+moinfo.sopi[h]; q++) {

        for(i=mo_offset[h]; i < mo_offset[h]+aoccpi[h]; i++)
          D_a[p][q] += C[p][i] * C[q][i];

        for(i=mo_offset[h]; i < mo_offset[h]+boccpi[h]; i++)
          D_b[p][q] += C[p][i] * C[q][i];
      }
    }
  }

  /* build the SO-basis Fock matrix */
  fock_a = block_matrix(nso, nso);
  fock_b = block_matrix(nso, nso);
  uhf_fock_build(fock_a, fock_b, D_a, D_b);
  free_block(D_a); free_block(D_b);

  /* transform the fock matrices to the MO bases */
  X = block_matrix(nso,nso);
  C_DGEMM('n','n',nso,nmo,nso,1.0,fock_a[0],nso,C[0],nmo,0.0,X[0],nso);
  C_DGEMM('t','n',nmo,nmo,nso,1.0,C[0],nmo,X[0],nso,0.0,fock_a[0],nso);

  C_DGEMM('n','n',nso,nmo,nso,1.0,fock_b[0],nso,C[0],nmo,0.0,X[0],nso);
  C_DGEMM('t','n',nmo,nmo,nso,1.0,C[0],nmo,X[0],nso,0.0,fock_b[0],nso);
  free_block(X);

  if(params.print_lvl > 2) {
    fprintf(outfile,"Alpha Fock matrix (before canonicalization)");
    print_mat(fock_a,nmo,nmo,outfile);
    fprintf(outfile,"Beta Fock matrix (before canonicalization)");
    print_mat(fock_b,nmo,nmo,outfile);
  }

  /** alpha Fock semicanonicalization **/

  Foo = (double ***) malloc(nirreps * sizeof(double **));
  Fvv = (double ***) malloc(nirreps * sizeof(double **));
  X = block_matrix(nmo, nmo);
  C_a = block_matrix(nso, nmo);
  alpha_evals = init_array(nmo);
  cnt = 0;
  for(h=0; h < nirreps; h++) {

    Foo[h] = block_matrix(aoccpi[h], aoccpi[h]);
    Fvv[h] = block_matrix(avirtpi[h], avirtpi[h]);

    for(i=mo_offset[h],I=0; i < mo_offset[h]+aoccpi[h]; i++,I++)
      for(j=mo_offset[h],J=0; j < mo_offset[h]+aoccpi[h]; j++,J++)
        Foo[h][I][J] = fock_a[i][j];

    for(a=mo_offset[h]+aoccpi[h],A=0; a < mo_offset[h]+moinfo.mopi[h]; a++,A++)
      for(b=mo_offset[h]+aoccpi[h],B=0; b < mo_offset[h]+moinfo.mopi[h]; b++,B++)
        Fvv[h][A][B] = fock_a[a][b];

    if(aoccpi[h]) {
      evals = init_array(aoccpi[h]);
      work = init_array(3*aoccpi[h]);
      if(stat = C_DSYEV('v','u', aoccpi[h], Foo[h][0], aoccpi[h], evals, work, aoccpi[h]*3)) {
        fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", h, stat);
        throw PsiException("transqt2: semicanonicalization error", __FILE__, __LINE__);
      }
      for(i=0; i<aoccpi[h]; i++) alpha_evals[cnt++] = evals[i];
      free(evals);
      free(work);

      for(i=mo_offset[h],I=0; i < mo_offset[h]+aoccpi[h]; i++,I++)
        for(j=mo_offset[h],J=0; j < mo_offset[h]+aoccpi[h]; j++,J++)
          X[i][j] = Foo[h][J][I];
    }

    if(avirtpi[h]) {
      evals = init_array(avirtpi[h]);
      work = init_array(3*avirtpi[h]);
      if(stat = C_DSYEV('v','u', avirtpi[h], Fvv[h][0], avirtpi[h], evals, work, avirtpi[h]*3)) {
        fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", h, stat);
        throw PsiException("transqt2: semicanonicalization error", __FILE__, __LINE__);
      }
      for(i=0; i<avirtpi[h]; i++) alpha_evals[cnt++] = evals[i];
      free(evals);
      free(work);

      for(a=mo_offset[h]+aoccpi[h],A=0; a < mo_offset[h]+moinfo.mopi[h]; a++,A++)
        for(b=mo_offset[h]+aoccpi[h],B=0; b < mo_offset[h]+moinfo.mopi[h]; b++,B++)
          X[a][b] = Fvv[h][B][A];
    }

    free_block(Foo[h]);
    free_block(Fvv[h]);
  }
  free(Foo);
  free(Fvv);
  free_block(fock_a);

  C_DGEMM('n','n',nso,nmo,nmo,1.0,C[0],nmo,X[0],nmo,0.0,C_a[0],nmo);

  free_block(X);

  /** beta Fock semicanonicalization **/

  Foo = (double ***) malloc(nirreps * sizeof(double **));
  Fvv = (double ***) malloc(nirreps * sizeof(double **));
  X = block_matrix(nmo, nmo);
  C_b = block_matrix(nso, nmo);
  beta_evals = init_array(nmo);
  cnt = 0;
  for(h=0; h < nirreps; h++) {
    /* leave the frozen core orbitals alone */
    for(i=mo_offset[h]; i < mo_offset[h]+moinfo.frdocc[h]; i++) X[i][i] = 1.0;

    Foo[h] = block_matrix(boccpi[h], boccpi[h]);
    Fvv[h] = block_matrix(bvirtpi[h], bvirtpi[h]);

    for(i=mo_offset[h],I=0; i < mo_offset[h]+boccpi[h]; i++,I++)
      for(j=mo_offset[h],J=0; j < mo_offset[h]+boccpi[h]; j++,J++)
        Foo[h][I][J] = fock_b[i][j];

    for(a=mo_offset[h]+boccpi[h],A=0; a < mo_offset[h]+moinfo.mopi[h]; a++,A++)
      for(b=mo_offset[h]+boccpi[h],B=0; b < mo_offset[h]+moinfo.mopi[h]; b++,B++)
        Fvv[h][A][B] = fock_b[a][b];

    if(boccpi[h]) {
      evals = init_array(boccpi[h]);
      work = init_array(3*boccpi[h]);
      if(stat = C_DSYEV('v','u', boccpi[h], Foo[h][0], boccpi[h], evals, work, boccpi[h]*3)) {
        fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", h, stat);
        throw PsiException("transqt2: semicanonicalization error", __FILE__, __LINE__);
      }
      for(i=0; i<boccpi[h]; i++) beta_evals[cnt++] = evals[i];
      free(evals);
      free(work);

      for(i=mo_offset[h],I=0; i < mo_offset[h]+boccpi[h]; i++,I++)
        for(j=mo_offset[h],J=0; j < mo_offset[h]+boccpi[h]; j++,J++)
          X[i][j] = Foo[h][J][I];
    }

    if(bvirtpi[h]) {
      evals = init_array(bvirtpi[h]);
      work = init_array(3*bvirtpi[h]);
      if(stat = C_DSYEV('v','u', bvirtpi[h], Fvv[h][0], bvirtpi[h], evals, work, bvirtpi[h]*3)) {
        fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", h, stat);
        throw PsiException("transqt2: semicanonicalization error", __FILE__, __LINE__);
      }
      for(i=0; i<bvirtpi[h]; i++) beta_evals[cnt++] = evals[i];
      free(evals);
      free(work);

      for(a=mo_offset[h]+boccpi[h],A=0; a < mo_offset[h]+moinfo.mopi[h]; a++,A++)
        for(b=mo_offset[h]+boccpi[h],B=0; b < mo_offset[h]+moinfo.mopi[h]; b++,B++)
          X[a][b] = Fvv[h][B][A];
    }

    free_block(Foo[h]);
    free_block(Fvv[h]);
  }
  free(Foo);
  free(Fvv);
  free_block(fock_b);

  C_DGEMM('n','n',nso,nmo,nmo,1.0,C[0],nmo,X[0],nmo,0.0,C_b[0],nmo);

  free_block(X);

  /* Write Semicanonical Alpha and Beta Fock Matrix Eigenvectors
     and Eigenvalues to the Checkpoint File */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_alpha_evals(alpha_evals);
  chkpt_wt_beta_evals(beta_evals);
  chkpt_wt_alpha_scf(C_a);
  chkpt_wt_beta_scf(C_b);
  chkpt_close();

  if(params.print_lvl > 2) {
    fprintf(outfile, "\nAlpha Eigenvalues\n");
    for (i=0; i<nmo; i++)
      fprintf(outfile, "%10.7lf\n", alpha_evals[i]);
    fprintf(outfile, "\nBeta Eigenvalues\n");
    for (i=0; i<nmo; i++)
      fprintf(outfile, "%10.7lf\n", beta_evals[i]);
    fflush(outfile);

    fprintf(outfile, "\nAlpha Eigenvectors\n");
    print_mat(C_a, nso, nmo, outfile);
    fprintf(outfile, "\nBeta Eigenvectors\n");
    print_mat(C_b, nso, nmo, outfile);
  }

  free_block(C_a);
  free_block(C_b);
  free_block(C);
  free(alpha_evals);
  free(beta_evals);
  free(mo_offset);
  free(so_offset);
}


void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b)
{
  int i, j, ij;
  int nso=0;
  int ntri=0;
  int stat=0;
  double *scratch;
  int lastbuf, idx;
  int p=0; int q=0; int r=0; int s=0;
  int pq=0; int rs=0;
  double value;
  Value *valptr;
  Label *lblptr;
  struct iwlbuf InBuf;
  double **Dt;

  nso = moinfo.nso;
  ntri = nso*(nso+1)/2;

  Dt = block_matrix(nso, nso);
  for(p=0; p < nso; p++)
    for(q=0; q < nso; q++)
      Dt[p][q] = D_a[p][q] + D_b[p][q];

  /* one-electron contributions */
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_T, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++) {
      fock_a[i][j] = fock_a[j][i] = scratch[ij];
      fock_b[i][j] = fock_b[j][i] = scratch[ij];
    }
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_V, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++) {
      fock_a[i][j] += scratch[ij];
      if(i!=j) fock_a[j][i] += scratch[ij];
      fock_b[i][j] += scratch[ij];
      if(i!=j) fock_b[j][i] += scratch[ij];
    }
  free(scratch);

  iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
  do {

    lastbuf = InBuf.lastbuf;
    lblptr = InBuf.labels;
    valptr = InBuf.values;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];
      value = (double) valptr[InBuf.idx];

      pq = INDEX(p,q);
      rs = INDEX(r,s);

      /* fprintf(outfile, "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value); */

      /* (pq|rs) */
      fock_a[p][q] += Dt[r][s] * value;
      fock_a[p][r] -= D_a[q][s] * value;
      fock_b[p][q] += Dt[r][s] * value;
      fock_b[p][r] -= D_b[q][s] * value;

      if(p!=q && r!=s && pq != rs) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (qp|sr) */
        fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
        fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (rs|qp) */
        fock_a[r][s] += Dt[q][p] * value;
        fock_a[r][q] -= D_a[s][p] * value;
        fock_b[r][s] += Dt[q][p] * value;
        fock_b[r][q] -= D_b[s][p] * value;

        /* (sr|pq) */
        fock_a[s][r] += Dt[p][q] * value;
        fock_a[s][p] -= D_a[r][q] * value;
        fock_b[s][r] += Dt[p][q] * value;
        fock_b[s][p] -= D_b[r][q] * value;

        /* (sr|qp) */
        fock_a[s][r] += Dt[q][p] * value;
        fock_a[s][q] -= D_a[r][p] * value;
        fock_b[s][r] += Dt[q][p] * value;
        fock_b[s][q] -= D_b[r][p] * value;
      }
      else if(p!=q && r!=s && pq==rs) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (qp|sr) */
        fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
        fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

      }
      else if(p!=q && r==s) {

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (rs|qp) */
        fock_a[r][s] += Dt[q][p] * value;
        fock_a[r][q] -= D_a[s][p] * value;
        fock_b[r][s] += Dt[q][p] * value;
        fock_b[r][q] -= D_b[s][p] * value;

      }
      else if(p==q && r!=s) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (sr|pq) */
        fock_a[s][r] += Dt[p][q] * value;
        fock_a[s][p] -= D_a[r][q] * value;
        fock_b[s][r] += Dt[p][q] * value;
        fock_b[s][p] -= D_b[r][q] * value;

      }
      else if(p==q && r==s && pq!=rs) {

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

      }
    }

    if(!lastbuf) iwl_buf_fetch(&InBuf);

  } while (!lastbuf);
  iwl_buf_close(&InBuf, 1);

  free_block(Dt);
}

  } // namespace transqt2
} // namespace psi
