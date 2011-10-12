/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"


namespace psi {
extern FILE *outfile;

namespace transqt {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

void semicanonical_fock(int averaged) /* averaged==0, regular semicanonical; =1 Z-averaged orbitals */
{
  int h, p, q, i, j, a, b, I, J, A, B, x, y, XX, YY;
  int stat, cnt;
  int *offset = init_int_array(moinfo.nirreps);
  int *aoccpi = init_int_array(moinfo.nirreps);
  int *asoccpi = init_int_array(moinfo.nirreps); /* socc per irrep, used for ZAPT */
  int *boccpi = init_int_array(moinfo.nirreps);
  int *avirtpi = init_int_array(moinfo.nirreps);
  int *bvirtpi = init_int_array(moinfo.nirreps);
  double **D_a = block_matrix(moinfo.nso, moinfo.nso);
  double **D_b = block_matrix(moinfo.nso, moinfo.nso);
  double **fock_a;
  double **fock_b;
  double **X;
  double ***Foo, ***Fvv, ***Fss; /* Fss for single-single block of Fock */
  double *evals, *work;
  double *alpha_evals;
  double *beta_evals;
  double **scf_vector_alpha;
  double **scf_vector_beta;

  for(h=0; h < moinfo.nirreps; h++) {
    if(averaged==1) { /* ZAPT */
        aoccpi[h] = moinfo.clsdpi[h];
        asoccpi[h] = moinfo.openpi[h];
    } else {
        aoccpi[h] = moinfo.clsdpi[h] + moinfo.openpi[h];
    }

    boccpi[h] = moinfo.clsdpi[h];
    avirtpi[h] = moinfo.virtpi[h];
    bvirtpi[h] = moinfo.virtpi[h] + moinfo.openpi[h];
  }

  for(h=1; h < moinfo.nirreps; h++) offset[h] = offset[h-1] + moinfo.orbspi[h-1];

  /* build the SO-basis alpha and beta densities from the ROHF eigenvector */
  for(h=0; h < moinfo.nirreps; h++) {
    for(p=offset[h]; p < offset[h]+moinfo.orbspi[h]; p++) {
      for(q=offset[h]; q < offset[h]+moinfo.orbspi[h]; q++) {
        if(averaged==1) {
          for(i=offset[h]; i < offset[h]+aoccpi[h]+asoccpi[h]; i++) {
            D_a[p][q] += moinfo.scf_vector[p][i] * moinfo.scf_vector[q][i];
          }
        } else {
          for(i=offset[h]; i < offset[h]+aoccpi[h]; i++) {
            D_a[p][q] += moinfo.scf_vector[p][i] * moinfo.scf_vector[q][i];
          }
        }
        for(i=offset[h]; i < offset[h]+boccpi[h]; i++) {
          D_b[p][q] += moinfo.scf_vector[p][i] * moinfo.scf_vector[q][i];
        }
      }
    }
  }

  /* build the SO-basis Fock matrix */
  fock_a = block_matrix(moinfo.nso, moinfo.nso);
  fock_b = block_matrix(moinfo.nso, moinfo.nso);
  uhf_fock_build(fock_a, fock_b, D_a, D_b);
  free_block(D_a); free_block(D_b);

  /* transform the fock matrices to the MO bases */
  X = block_matrix(moinfo.nso,moinfo.nso);
  C_DGEMM('n','n',moinfo.nso,moinfo.nmo,moinfo.nso,1.0,&(fock_a[0][0]),moinfo.nso,&(moinfo.scf_vector[0][0]),moinfo.nmo,
          0,&(X[0][0]),moinfo.nso);
  C_DGEMM('t','n',moinfo.nmo,moinfo.nmo,moinfo.nso,1.0,&(moinfo.scf_vector[0][0]),moinfo.nmo,&(X[0][0]),moinfo.nso,
          0,&(fock_a[0][0]),moinfo.nso);

  C_DGEMM('n','n',moinfo.nso,moinfo.nmo,moinfo.nso,1.0,&(fock_b[0][0]),moinfo.nso,&(moinfo.scf_vector[0][0]),moinfo.nmo,
          0,&(X[0][0]),moinfo.nso);
  C_DGEMM('t','n',moinfo.nmo,moinfo.nmo,moinfo.nso,1.0,&(moinfo.scf_vector[0][0]),moinfo.nmo,&(X[0][0]),moinfo.nso,
          0,&(fock_b[0][0]),moinfo.nso);
  free_block(X);


  /** alpha Fock semicanonicalization **/

  Foo = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  Fvv = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  if(averaged==1)
    Fss = (double ***) malloc(moinfo.nirreps * sizeof(double **));

  X = block_matrix(moinfo.nmo, moinfo.nmo);
  scf_vector_alpha = block_matrix(moinfo.nmo, moinfo.nmo);
  alpha_evals = init_array(moinfo.nmo);
  cnt = 0;
  for(h=0; h < moinfo.nirreps; h++) {

    Foo[h] = block_matrix(aoccpi[h], aoccpi[h]);
    Fvv[h] = block_matrix(avirtpi[h], avirtpi[h]);
   if(averaged==1)
        Fss[h] = block_matrix(asoccpi[h], asoccpi[h]);


    if(averaged==1) { 		/* Average Alpha and Beta Fock Matrices */
      for(i=offset[h],I=0; i < offset[h]+aoccpi[h]; i++,I++)
        for(j=offset[h],J=0; j < offset[h]+aoccpi[h]; j++,J++)
          Foo[h][I][J] = (fock_a[i][j] + fock_b[i][j])/2.0;

      for(x=offset[h]+aoccpi[h],XX=0; x < offset[h]+aoccpi[h]+asoccpi[h]; x++,XX++)
        for(y=offset[h]+aoccpi[h],YY=0; y < offset[h]+aoccpi[h]+asoccpi[h]; y++,YY++)
          Fss[h][XX][YY] = (fock_a[x][y] + fock_b[x][y])/2.0;

      for(a=offset[h]+aoccpi[h]+asoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
        for(b=offset[h]+aoccpi[h]+asoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
          Fvv[h][A][B] = (fock_a[a][b] + fock_b[a][b])/2.0;
    } else {
      for(i=offset[h],I=0; i < offset[h]+aoccpi[h]; i++,I++)
        for(j=offset[h],J=0; j < offset[h]+aoccpi[h]; j++,J++)
          Foo[h][I][J] = fock_a[i][j];

      for(a=offset[h]+aoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
        for(b=offset[h]+aoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
          Fvv[h][A][B] = fock_a[a][b];
    }

    /* diagonalize in docc-docc subspace */
    if(aoccpi[h]) {
      evals = init_array(aoccpi[h]);
      work = init_array(3*aoccpi[h]);
      if(stat = C_DSYEV('v','u', aoccpi[h], &(Foo[h][0][0]), aoccpi[h], evals, work, aoccpi[h]*3)) {

        fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n",
                h, stat);
        exit(PSI_RETURN_FAILURE);
      }
      for(i=0; i<aoccpi[h]; i++) alpha_evals[cnt++] = evals[i];
      free(evals);
      free(work);

    for(i=offset[h],I=0; i < offset[h]+aoccpi[h]; i++,I++)
      for(j=offset[h],J=0; j < offset[h]+aoccpi[h]; j++,J++)
        X[i][j] = Foo[h][J][I];
    }

    /* diagonalize in socc-socc subspace */
    if(averaged==1 && asoccpi[h]) {
      evals = init_array(asoccpi[h]);
      work = init_array(3*asoccpi[h]);
      if(stat = C_DSYEV('v','u', asoccpi[h], &(Fss[h][0][0]), asoccpi[h], evals, work, asoccpi[h]*3)) {
        fprintf(outfile, "rotate(): Error in alpha Fss[%1d] diagonalization. stat = %d\n", h, stat);
        exit(PSI_RETURN_FAILURE);
      }
      for(i=0; i<asoccpi[h]; i++) alpha_evals[cnt++] = evals[i];
      free(evals);
      free(work);

      for(x=offset[h]+aoccpi[h],XX=0; x < offset[h]+aoccpi[h]+asoccpi[h]; x++,XX++)
        for(y=offset[h]+aoccpi[h],YY=0; y < offset[h]+aoccpi[h]+asoccpi[h]; y++,YY++)
          X[x][y] = Fss[h][YY][XX];
    }


    /* diagonalize in vir-vir subspace */
    if(avirtpi[h]) {
      evals = init_array(avirtpi[h]);
      work = init_array(3*avirtpi[h]);
      if(stat = C_DSYEV('v','u', avirtpi[h], &(Fvv[h][0][0]), avirtpi[h],
                        evals, work, avirtpi[h]*3)) {
      fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n",
              h, stat);
      exit(PSI_RETURN_FAILURE);
    }
    for(i=0; i<avirtpi[h]; i++) alpha_evals[cnt++] = evals[i];
    free(evals);
    free(work);

      if(averaged==1) {
        for(a=offset[h]+aoccpi[h]+asoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
          for(b=offset[h]+aoccpi[h]+asoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
            X[a][b] = Fvv[h][B][A];
      } else {
        for(a=offset[h]+aoccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
          for(b=offset[h]+aoccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
            X[a][b] = Fvv[h][B][A];
      }

    }

    free_block(Foo[h]);
    free_block(Fvv[h]);
    if(averaged==1)
        free_block(Fss[h]);
  }
  free(Foo);
  free(Fvv);
  if(averaged==1)
    free(Fss);
  free_block(fock_a);

  C_DGEMM('n','n',moinfo.nso,moinfo.nmo,moinfo.nmo,1,&(moinfo.scf_vector[0][0]),moinfo.nmo,&(X[0][0]),
          moinfo.nmo,0,&(scf_vector_alpha[0][0]),moinfo.nmo);

  free_block(X);

  /** beta Fock semicanonicalization **/
  if(averaged==0) { /* No beta Fock matrix for ZAPT */

    Foo = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    Fvv = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    X = block_matrix(moinfo.nmo, moinfo.nmo);
    scf_vector_beta = block_matrix(moinfo.nmo, moinfo.nmo);
    beta_evals = init_array(moinfo.nmo);
    cnt = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      /* leave the frozen core orbitals alone */
      for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++) X[i][i] = 1.0;

      Foo[h] = block_matrix(boccpi[h], boccpi[h]);
      Fvv[h] = block_matrix(bvirtpi[h], bvirtpi[h]);

      for(i=offset[h],I=0; i < offset[h]+boccpi[h]; i++,I++)
        for(j=offset[h],J=0; j < offset[h]+boccpi[h]; j++,J++)
          Foo[h][I][J] = fock_b[i][j];

        for(a=offset[h]+boccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
          for(b=offset[h]+boccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
            Fvv[h][A][B] = fock_b[a][b];

        if(boccpi[h]) {
          evals = init_array(boccpi[h]);
          work = init_array(3*boccpi[h]);
          if(stat = C_DSYEV('v','u', boccpi[h], &(Foo[h][0][0]),
                          boccpi[h], evals, work, boccpi[h]*3)) {
            fprintf(outfile, "rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n",
                  h, stat);
            exit(PSI_RETURN_FAILURE);
          }
          for(i=0; i<boccpi[h]; i++) beta_evals[cnt++] = evals[i];
          free(evals);
          free(work);

          for(i=offset[h],I=0; i < offset[h]+boccpi[h]; i++,I++)
            for(j=offset[h],J=0; j < offset[h]+boccpi[h]; j++,J++)
              X[i][j] = Foo[h][J][I];
        }

        if(bvirtpi[h]) {
          evals = init_array(bvirtpi[h]);
          work = init_array(3*bvirtpi[h]);
          if(stat = C_DSYEV('v','u', bvirtpi[h], &(Fvv[h][0][0]), bvirtpi[h],
                          evals, work, bvirtpi[h]*3)) {
            fprintf(outfile, "rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n",
                  h, stat);
            exit(PSI_RETURN_FAILURE);
          }
          for(i=0; i<bvirtpi[h]; i++) beta_evals[cnt++] = evals[i];
          free(evals);
          free(work);

          for(a=offset[h]+boccpi[h],A=0; a < offset[h]+moinfo.orbspi[h]; a++,A++)
            for(b=offset[h]+boccpi[h],B=0; b < offset[h]+moinfo.orbspi[h]; b++,B++)
              X[a][b] = Fvv[h][B][A];
       }

        free_block(Foo[h]);
        free_block(Fvv[h]);
      }
      free(Foo);
      free(Fvv);
      free_block(fock_b);

    C_DGEMM('n','n',moinfo.nso,moinfo.nmo,moinfo.nmo,1,&(moinfo.scf_vector[0][0]),moinfo.nmo,&(X[0][0]),
          moinfo.nmo,0,&(scf_vector_beta[0][0]),moinfo.nmo);

    free_block(X);
  }

  /* Write Semicanonical Alpha and Beta Fock Matrix Eigenvectors
     and Eigenvalues to the Checkpoint File */
  if(averaged == 0) {
    chkpt_wt_alpha_evals(alpha_evals);
    chkpt_wt_alpha_scf(scf_vector_alpha);
    chkpt_wt_beta_evals(beta_evals);
    chkpt_wt_beta_scf(scf_vector_beta);
  } else {
    chkpt_wt_evals(alpha_evals);
    chkpt_wt_scf(scf_vector_alpha);
  }


  /*fprintf(outfile,"\nAlpha Eigenvalues\n");
  for(i=0; i<moinfo.nmo; i++) fprintf(outfile,"%10.7lf\n",alpha_evals[i]);
  fprintf(outfile,"\nBeta Eigenvalues\n");
  for(i=0; i<moinfo.nmo; i++) fprintf(outfile,"%10.7lf\n",beta_evals[i]);
  fflush(outfile);

  fprintf(outfile,"\nAlpha Eigenvalues\n");
  print_mat(scf_vector_alpha,moinfo.nmo,moinfo.nmo,outfile);
  fprintf(outfile,"\nBeta Eigenvalues\n");
  print_mat(scf_vector_beta,moinfo.nmo,moinfo.nmo,outfile);*/

  free_block(scf_vector_alpha);
  free(alpha_evals);
  if(averaged==0) {
    free_block(scf_vector_beta);
    free(beta_evals);
  }
  free(offset);
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
  scratch = init_array(moinfo.noeints);
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

}} // end namespace psi::transqt
