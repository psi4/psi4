/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

/* the non-symmetric transpose cases have not been tested */

int dpd_dot13(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
    int transt, int transz, double alpha, double beta)
{
  int h, Gp, Gq, Gr, Gs, GT, GI, GZ, Tblock, Zblock;
  int p, q, r, s;
  int P, Q, R, S;
  int row, col;
  int nirreps;
  double **X;
  double value;

  nirreps = T->params->nirreps;
  GT = T->my_irrep;
  GI = I->file.my_irrep;
  GZ = Z->my_irrep;

  /* Get the two-index quantities from disk */
  dpd_file2_mat_init(T);
  dpd_file2_mat_rd(T);

  dpd_file2_scm(Z, beta);
  dpd_file2_mat_init(Z);
  dpd_file2_mat_rd(Z);

#ifdef DPD_TIMER
  timer_on("dot13");
#endif

  /* loop over row irreps of bufI; h = Gpq = Grs^GI */
  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(I, h);
    dpd_buf4_mat_irrep_rd(I, h);

    /* Loop over row irreps of target Z, Gz = Gqs */
    for(Gq=0; Gq < nirreps; Gq++) {
      /* Gs = Gq;  Gp = Gr = h^Gq; */
      Gp = h^Gq; Gr = h^Gq^GT; Gs = Gq^GZ;
      if (!transt) Tblock = Gp; else Tblock = Gr;
      if (!transz) Zblock = Gq; else Zblock = Gs;

      /* Allocate space for the X buffer */
      if(T->params->ppi[Gp] && T->params->qpi[Gr])
        X = dpd_block_matrix(T->params->ppi[Gp],T->params->qpi[Gr]);

      /* Loop over orbitals of the target */
      for(q=0; q < Z->params->ppi[Gq]; q++) {
        Q = Z->params->poff[Gq] + q;
        for(s=0; s < Z->params->qpi[Gs]; s++) {
          S = Z->params->qoff[Gs] + s;

          /* Loop over orbitals of the two-index term */
          for(p=0; p < T->params->ppi[Gp]; p++) {
            P = T->params->poff[Gp] + p;
            for(r=0; r < T->params->qpi[Gr]; r++) {
              R = T->params->qoff[Gr] + r;

              /* Calculate row and column indices in I */
              if(!transt && !transz) {
                row = I->params->rowidx[P][Q];
                col = I->params->colidx[R][S];
              }
              else if(transt && !transz) {
                row = I->params->rowidx[R][Q];
                col = I->params->colidx[P][S];
              }
              else if(!transt && transz) {
                row = I->params->rowidx[P][S];
                col = I->params->colidx[R][Q];
              }
              else if(transt && transz) {
                row = I->params->rowidx[R][S];
                col = I->params->colidx[P][Q];
              }

              /* Build the X buffer */
              X[p][r] = I->matrix[h][row][col]; 

            }
          }

          value = dot_block(T->matrix[Tblock], X, T->params->ppi[Gp],
              T->params->qpi[Gr], alpha); 

          Z->matrix[Zblock][q][s] += value;
        }
      }
      if(T->params->ppi[Gp] && T->params->qpi[Gr])
        dpd_free_block(X, T->params->ppi[Gp],T->params->qpi[Gr]);
    }
    dpd_buf4_mat_irrep_close(I, h);
  }

#ifdef DPD_TIMER
  timer_off("dot13");
#endif

  /* Close the two-index quantities */
  dpd_file2_mat_close(T);
  dpd_file2_mat_wrt(Z);
  dpd_file2_mat_close(Z);

  return 0;
}

} // namespace psi
