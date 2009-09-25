/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

#define HALFPI 1.5707963268

/*
** Follow the eigenvector of the MO Hessian downhill.  This version
** is based on orbital rotations as found in the detcas code.  It 
** appears to be working better than the antisymmetric matrix
** alternative found below.
**
** CDS 8/23/03
*/
void follow_evec_UHF(double *vf, int dim_A, int dim_B)
{
  int dim, h, nirreps;
  int i, ii, a, aa;
  int nso, nmo;
  double *v_A, *v_B, theta, sintheta, costheta, **scf, scale;

  chkpt_init(PSIO_OPEN_OLD);

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  nirreps = moinfo.nirreps;

  dim = dim_A + dim_B;
  v_A = vf;
  v_B = vf + dim_A;
 
  scale = HALFPI;  /* full step would be Pi/2 */
  scale *= params.scale;  
 
  /* rotate the alpha orbitals */
  scf = chkpt_rd_alpha_scf();

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tOld molecular orbitals (alpha):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  for (h=0; h<nirreps; h++) {
    for (a=0; a<moinfo.avirtpi[h]; a++) {
      /* to rotate the C matrix we need to convert indices to Piter order */
      aa = moinfo.qt2pitzer_a[moinfo.qt_avir[a] + moinfo.avir_off[h]];
      /* i must have same irrep if rotation is totally symmetric */ 
      for (i=0; i<moinfo.aoccpi[h]; i++) {
        ii = moinfo.qt2pitzer_a[moinfo.qt_aocc[i] + moinfo.aocc_off[h]];
        /* The sign on theta doesn't appear to matter, which is good... */
        theta = *v_A++ * scale;
        costheta = cos(theta);
        sintheta = sin(theta);
        C_DROT(nso,&(scf[0][ii]),nso,&(scf[0][aa]),nso,costheta,sintheta);
      } /* end loop over i */
    } /* end loop over a */
  }  /* end loop over h */

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tNew molecular orbitals (alpha):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  chkpt_wt_alpha_scf(scf);
  free_block(scf);


  /* rotate the beta orbitals */
  scf = chkpt_rd_beta_scf();

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tOld molecular orbitals (beta):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  for (h=0; h<nirreps; h++) {
    for (a=0; a<moinfo.bvirtpi[h]; a++) {
      /* to rotate the C matrix we need to convert indices to Piter order */
      aa = moinfo.qt2pitzer_b[moinfo.qt_bvir[a] + moinfo.bvir_off[h]];
      /* i must have same irrep if rotation is totally symmetric */ 
      for (i=0; i<moinfo.boccpi[h]; i++) {
        ii = moinfo.qt2pitzer_b[moinfo.qt_bocc[i] + moinfo.bocc_off[h]];
        theta = *v_B++ * scale;
        /* The sign on theta doesn't appear to matter, which is good... */
        costheta = cos(theta);
        sintheta = sin(theta);
        C_DROT(nso,&(scf[0][ii]),nso,&(scf[0][aa]),nso,costheta,sintheta);
      } /* end loop over i */
    } /* end loop over a */
  }  /* end loop over h */

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tNew molecular orbitals (beta):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  chkpt_wt_beta_scf(scf);
  free_block(scf);

  chkpt_close();

}


/* 
** Just for fun, compare this alternative way to get the rotation.
** This version is based on Daniel's orbital rotation code in
** ccenergy/rotate.c
**
** CDS 8/23/03
*/
void follow_evec_UHF2(double *vf, int dim_A, int dim_B)
{
  int dim, h, nirreps;
  int i, ii, a, aa, ai, j, ij;
  int nso, nmo, ntri, stat;
  double *v_A, *v_B, **scf, **scf_new;
  double **scf_a_orig, **scf_b_orig, **scf_a, **scf_b;
  double **U, **MO_S, **SO_S, **X, *evals, *work, **S;
  double *scratch;

  chkpt_init(PSIO_OPEN_OLD);

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  nirreps = moinfo.nirreps;

  dim = dim_A + dim_B;
  v_A = vf;
  v_B = vf + dim_A;

  SO_S = block_matrix(nso, nso);
  ntri = nso * (nso+1)/2;
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_S, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      SO_S[i][j] = SO_S[j][i] = scratch[ij];
    }
  free(scratch);
 
  /* rotate the alpha orbitals */
  U = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) U[i][i] = 1.0;

  scf = chkpt_rd_alpha_scf();
  scf_a_orig = chkpt_rd_alpha_scf();
  scf_new = block_matrix(nso, nmo);

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tOld molecular orbitals (alpha):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  for (h=0,ai=0; h<nirreps; h++) {
    for (a=0; a<moinfo.avirtpi[h]; a++) {
      /* to rotate the C matrix we need to convert indices to Piter order */
      aa = moinfo.qt2pitzer_a[moinfo.qt_avir[a] + moinfo.avir_off[h]];
      /* i must have same irrep if rotation is totally symmetric */ 
      for (i=0; i<moinfo.aoccpi[h]; i++,ai++) {
        ii = moinfo.qt2pitzer_a[moinfo.qt_aocc[i] + moinfo.aocc_off[h]];
        U[ii][aa] = v_A[ai];
        U[aa][ii] = -v_A[ai];
      } /* end loop over i */
    } /* end loop over a */
  }  /* end loop over h */

  C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
            0,&(scf_new[0][0]),nmo);
  free_block(U);
  free_block(scf);

  MO_S = block_matrix(nmo, nmo);

  /* transform the overlap into the new alpha MO basis */
  X = block_matrix(nso, nso);
  C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]),
          nso, 0, &(X[0][0]), nso);
  C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]),
          nmo, 0, &(MO_S[0][0]), nmo);
  free_block(X);
                                                                               
  evals = init_array(nmo);
  work = init_array(nmo*3);
  if(stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3)) {
    fprintf(outfile, 
      "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
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

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tNew molecular orbitals (alpha):\n");
    print_mat(scf_a, nso, nmo, outfile);
  }

  chkpt_wt_alpha_scf(scf_a);
  free_block(scf_a);


  /* rotate the beta orbitals */
  U = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) U[i][i] = 1.0;

  scf = chkpt_rd_beta_scf();
  scf_b_orig = chkpt_rd_beta_scf();

  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tOld molecular orbitals (beta):\n");
    print_mat(scf, nso, nmo, outfile);
  }

  for (h=0,ai=0; h<nirreps; h++) {
    for (a=0; a<moinfo.bvirtpi[h]; a++) {
      /* to rotate the C matrix we need to convert indices to Piter order */
      aa = moinfo.qt2pitzer_b[moinfo.qt_bvir[a] + moinfo.bvir_off[h]];
      /* i must have same irrep if rotation is totally symmetric */ 
      for (i=0; i<moinfo.boccpi[h]; i++,ai++) {
        ii = moinfo.qt2pitzer_b[moinfo.qt_bocc[i] + moinfo.bocc_off[h]];
        U[ii][aa] = v_B[ai];
        U[aa][ii] = -v_B[ai];        
      } /* end loop over i */
    } /* end loop over a */
  }  /* end loop over h */

  scf_new = block_matrix(nso, nmo);
  C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
          0,&(scf_new[0][0]),nmo);
  free_block(U);
  free_block(scf);
                                                                                
  MO_S = block_matrix(nmo, nmo);
                                                                                
  /* transform the overlap into the new beta MO basis */
  X = block_matrix(nso, nso);
  C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]),
          nso, 0, &(X[0][0]), nso);
  C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]),
          nmo, 0, &(MO_S[0][0]), nmo);
  free_block(X);
                                                                                
  evals = init_array(nmo);
  work = init_array(nmo*3);
  if(stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3)) {
    fprintf(outfile, 
      "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
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


  if (params.print_lvl > 2) {
    fprintf(outfile, "\n\tNew molecular orbitals (beta):\n");
    print_mat(scf_b, nso, nmo, outfile);
  }

  chkpt_wt_beta_scf(scf_b);
  free_block(scf_b);

  free_block(SO_S);
  chkpt_close();

}



}} // namespace psi::stable
