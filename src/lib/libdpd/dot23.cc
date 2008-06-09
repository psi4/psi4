/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {
	
/* dpd_dot23(): Contracts a two-index-electron dpd with a four-index dpd
 ** where both indices in the former and indices two and three in the
 ** latter are summed.
 **
 ** Arguments:
 **   dpdfile2 *T: A pointer to the two-index-electron factor.
 **   dpdbuf4 *I: A pointer to the four-index factor.
 **   dpdfile2 *Z: A pointer to the two-index target.
 **   int transt: A boolean indicating whether the T-factor should be
 **               transposed.
 **   int transz: A boolean indicating whether the Z-product should be
 **               transposed.
 **   double alpha: A prefactor for the product alpha * T * I.
 **   double beta: A prefactor for the target beta * Z.
 **
 ** Example contractions:
 **    beta * Z(p,s) = alpha * T(q,r) * I(pq,rs)
 **           (transt = 0; transz =0;)
 **    beta * Z(p,s) = alpha * T(r,q) * I(pq,rs)
 **           (transt = 1; transz =0;)
 **    beta * Z(s,p) = alpha * T(q,r) * I(pq,rs)
 **           (transt = 0; transz =1;)
 **    beta * Z(s,p) = alpha * T(r,q) * I(pq,rs)
 **           (transt = 1; transz =1;)
 ** non-symmetric states have not been tested for transposed cases
 */


int dpd_dot23(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
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
  timer_on("dot23");
#endif

  /* loop over row irreps of buf I(pq,rs); h = Gpq = Grs^GI */
  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(I, h);
    dpd_buf4_mat_irrep_rd(I, h);

    /* Loop over irreps of the target Z; GZ = Gps */
    for(Gp=0; Gp < nirreps; Gp++) {
      /* Gs = Gp;  Gq = Gr = h^Gp; */
      Gq = h^Gp; Gr = h^Gp^GT; Gs = Gp^GZ;
      if (!transt) Tblock = Gq; else Tblock = Gr;
      if (!transz) Zblock = Gp; else Zblock = Gs;

      /* Allocate space for the X buffer */
      if(T->params->ppi[Gq] && T->params->qpi[Gr])
        X = dpd_block_matrix(T->params->ppi[Gq],T->params->qpi[Gr]);

      /* Loop over orbitals of the target */
      for(p=0; p < Z->params->ppi[Gp]; p++) {
        P = Z->params->poff[Gp] + p;
        for(s=0; s < Z->params->qpi[Gs]; s++) {
          S = Z->params->qoff[Gs] + s;

          /* Loop over orbitals of the two-index term */
          for(q=0; q < T->params->ppi[Gq]; q++) {
            Q = T->params->poff[Gq] + q;
            for(r=0; r < T->params->qpi[Gr]; r++) {
              R = T->params->qoff[Gr] + r;

              /* Calculate row and column indices in I */
              if(!transt && !transz) {
                row = I->params->rowidx[P][Q];
                col = I->params->colidx[R][S];
              }
              else if(transt && !transz) {
                row = I->params->rowidx[P][R];
                col = I->params->colidx[Q][S];
              }
              else if(!transt && transz) {
                row = I->params->rowidx[S][Q];
                col = I->params->colidx[R][P];
              }
              else if(transt && transz) {
                row = I->params->rowidx[S][R];
                col = I->params->colidx[Q][P];
              }

              /* Build the X buffer */
              X[q][r] = I->matrix[h][row][col]; 

            }
          }

          value = dot_block(T->matrix[Tblock], X, T->params->ppi[Gq],
              T->params->qpi[Gr], alpha); 

          Z->matrix[Zblock][p][s] += value;
        }
      }
      if(T->params->ppi[Gq] && T->params->qpi[Gr])
        dpd_free_block(X, T->params->ppi[Gq],T->params->qpi[Gr]);
    }
    dpd_buf4_mat_irrep_close(I, h);
  }

#ifdef DPD_TIMER
  timer_off("dot23");
#endif

  /* Close the two-index quantities */
  dpd_file2_mat_close(T);
  dpd_file2_mat_wrt(Z);
  dpd_file2_mat_close(Z);

  return 0;
}

} // namespace psi
