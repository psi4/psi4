/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_W(void);
void uhf_W(void);
void rhf_sf_relax_I(void);
void uhf_sf_relax_I(void);

void relax_I(void)
{
  if(params.ref == 0) rhf_sf_relax_I();
  else if(params.ref == 2) uhf_sf_relax_I();
}

void rhf_W(void)
{
  int i, j, a, b, h;
  int nirreps;
  double fii = 0.0;
  double fjj = 0.0;
  double faa = 0.0;
  double fbb = 0.0;
  double Dij = 0.0;
  double Dab = 0.0;
  double Dai = 0.0;
  dpdfile2 W;
  dpdfile2 D;
  dpdfile2 f;
  dpdbuf4 T;
  dpdbuf4 I;

  nirreps = mo.nirreps;

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 0, "WIJ");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract442(&T, &I, &W, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 1, "WAB");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract442(&T, &I, &W, 3, 3, -2.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 0, "WAI");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
  global_dpd_->contract442(&T, &I, &W, 2, 0, -4.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 0, "WIJ");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(i=0; i < mo.occpi[h]; i++) {
      fii = f.matrix[h][i][i];
      for(j=0; j < mo.occpi[h]; j++) {
        fjj = f.matrix[h][j][j];
        Dij = D.matrix[h][i][j];
        W.matrix[h][i][j] -= 0.5 * Dij * ( fii + fjj );
      }
    }
  }
  global_dpd_->file2_mat_wrt(&W);
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 1, "WAB");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      faa = f.matrix[h][a][a];
      for(b=0; b < mo.virpi[h]; b++) {
        fbb = f.matrix[h][b][b];
        Dab = D.matrix[h][a][b];
        W.matrix[h][a][b] -= 0.5 * Dab * ( faa + fbb );
      }
    }
  }

  global_dpd_->file2_mat_wrt(&W);
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 0, "WAI");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      for(i=0; i < mo.occpi[h]; i++) {
        fii = f.matrix[h][i][i];
        Dai = D.matrix[h][a][i];
        W.matrix[h][a][i] -= Dai * fii;
      }
    }
  }

  global_dpd_->file2_mat_wrt(&W);
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 1, "WIA");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      for(i=0; i < mo.occpi[h]; i++) {
        faa = f.matrix[h][a][a];
        Dai = D.matrix[h][a][i];
        W.matrix[h][i][a] += Dai * faa;
      }
    }
  }

  global_dpd_->file2_mat_wrt(&W);
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 0, "WIJ");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  global_dpd_->buf4_init(&I, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  global_dpd_->dot13(&D, &I, &W, 0, 0, -2.0, 1.0);
  global_dpd_->dot14(&D, &I, &W, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  global_dpd_->buf4_init(&I, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  global_dpd_->dot24(&D, &I, &W, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
  global_dpd_->dot23(&D, &I, &W, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&D);
  global_dpd_->buf4_close(&I);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "DAI");
  global_dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->dot13(&D, &I, &W, 0, 0, -2.0, 1.0);
  global_dpd_->dot13(&D, &I, &W, 0, 1, -2.0, 1.0);
  global_dpd_->dot14(&D, &I, &W, 0, 0, 1.0, 1.0);
  global_dpd_->dot14(&D, &I, &W, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&f);
  global_dpd_->file2_mat_rd(&f);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 0, "WIJ");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(i=0; i < mo.occpi[h]; i++) {
      fii = f.matrix[h][i][i];
      W.matrix[h][i][i] += 2.0 * fii;
    }
  }
  global_dpd_->file2_mat_wrt(&W);
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);
  global_dpd_->file2_mat_close(&f);
  global_dpd_->file2_close(&f);
}

void uhf_W(void)
{

}

void rhf_sf_relax_I(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, j, e, *occpi, *virtpi;

  nirreps = mo.nirreps;
  occpi = mo.occpi;
  virtpi = mo.virtpi;

  /*** occupied-virtual relaxation terms */

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,A)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /* I(i,a) = I'(i,a) - sum_m f(i,m) D(orb)(a,m) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,a)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(i,a)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_E,M D(orb)(E,M) [<EI||MJ> + <EJ||MI>]
                      - 2 sum_e,m D(orb)(e,m) <eI|mJ> */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,J)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);

  /* I(i,j) <-- I'(i,j) - sum_e,m D(orb)(e,m) [<ei||mj> + <ej||mi>]
                      - 2 sum_E,M D(orb)(E,M) <Ei|Mj> */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(i,j)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->dot13(&D, &E, &I, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -2.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);
}

void uhf_sf_relax_I(void)
{

}

}} /* End namespaces */
