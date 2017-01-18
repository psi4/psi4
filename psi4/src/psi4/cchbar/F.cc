/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/

/*! \defgroup CCHBAR cchbar: Compute the similarity-transformed Hamiltonian */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* F_build(): Constructs the one-electron HBAR matrix elements
** Fme, Fae, and Fmi.  These are defined in spin-orbitals as
**
** Fme = f_me + t_n^f <mn||ef>
**
** Fae = f_ae - 1/2 f_me t_m^a + f_m^f <am||ef> - 1/2 taut_mn^af <mn||ef>
**
** Fmi = f_mi + 1/2 f_me t_i^e + t_n^e <mn||ie> + 1/2 taut_in^ef <mn||ef>
**
** where taut_ij^ab = t_ij^ab + 1/2 ( t_i^a t_j^b - t_i^b t_j^a )
**
** The standard named FAE, Fae, FMI, and Fmi are used for the complete
** matrix elements, while the names FAEt, Faet, FMIt, and Fmit are used for
** matrices with the diagonal elements removed.
**
** TDC, revised June 2002
*/

void F_build(void) {
  int h,e,a,m;
  dpdfile2 Faet, FAEt, Fmit, FMIt;
  dpdfile2 Fae, FAE, FMI, Fmi;
  dpdfile2 fab, fAB, fij, fIJ;
  dpdfile2 FME, Fme;
  dpdfile2 fIA, fia;
  dpdfile2 tIA, tia;
  dpdbuf4 F_anti, F, E_anti, E, D_anti, D;
  dpdbuf4 tautIJAB, tautijab, tautIjAb, taut;
  int Gma, Gm, Ga, Gf, Ge, ma, M, A, fe, ef, f, nrows, ncols;
  double *X;

  if(params.ref == 0) {

    /** FME **/
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);

    /** FAE **/
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");

    /*
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F 2<ia|bc> - <ia|cb>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&F);
    */

    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&tIA);
    global_dpd_->file2_mat_rd(&tIA);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    for(Gma=0; Gma < moinfo.nirreps; Gma++) {
      global_dpd_->buf4_mat_irrep_row_init(&F, Gma);
      X = init_array(F.params->coltot[Gma]);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);
	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = Ge = Gm ^ Gma;  /* Fae is totally symmetric */
	Gf = Gm; /* T1 is totally symmetric */
	M = m - F.params->poff[Gm];
	A = a - F.params->qoff[Ga];

	zero_arr(X, F.params->coltot[Gma]);

	/* build spin-adapted F-integrals for current ma */
	for(fe=0; fe < F.params->coltot[Gma]; fe++) {
	  f = F.params->colorb[Gma][fe][0];
	  e = F.params->colorb[Gma][fe][1];
	  ef = F.params->colidx[e][f];
	  X[fe] = 2.0 * F.matrix[Gma][0][fe] - F.matrix[Gma][0][ef];
	}

	nrows = moinfo.virtpi[Gf];
	ncols = moinfo.virtpi[Ge];
	if(nrows && ncols)
	  C_DGEMV('t',nrows,ncols,1.0,&X[F.col_offset[Gma][Gf]],ncols,
		  tIA.matrix[Gm][M],1,1.0,
		  FAE.matrix[Ga][A],1);
      }

      free(X);
      global_dpd_->buf4_mat_irrep_row_close(&F, Gma);
    }
    global_dpd_->buf4_close(&F);
    global_dpd_->file2_mat_close(&tIA);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_mat_wrt(&FAE);
    global_dpd_->file2_mat_close(&FAE);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    global_dpd_->contract442(&tautIjAb, &D, &FAE, 3, 3, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tautIjAb);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->dot13(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_mat_init(&FAEt);
    global_dpd_->file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FAEt);
    global_dpd_->file2_mat_close(&FAEt);
    global_dpd_->file2_close(&FAEt);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_mat_init(&FMIt);
    global_dpd_->file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FMIt);
    global_dpd_->file2_mat_close(&FMIt);
    global_dpd_->file2_close(&FMIt);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* FME and Fme */
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    global_dpd_->dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

    /* FAE and Fae */

    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    global_dpd_->file2_close(&fab);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    global_dpd_->dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&F_anti);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    global_dpd_->contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    global_dpd_->contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&tautijab);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    global_dpd_->contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    global_dpd_->contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tautIjAb);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->contract222(&tia, &Fme, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&Fme);

    /* Form Fae and FAE tilde intermediates */
    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);


    /* FMI and Fmi */

    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    global_dpd_->file2_close(&fij);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    global_dpd_->dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0);

    global_dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&E_anti);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautijab");

    global_dpd_->contract442(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1);
    global_dpd_->contract442(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1);

    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&tautijab);
    global_dpd_->buf4_close(&D_anti);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->contract442(&D, &tautIjAb, &Fmi, 1, 1, 1, 1);

    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
    global_dpd_->contract222(&Fme, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&Fme);
    global_dpd_->file2_close(&tia);

    /* FMI and Fmi tilde intermediates */
    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);


    /* remove diagonal elements from Ft's */
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 1, 1, "Faet");
    global_dpd_->file2_mat_init(&Faet);
    global_dpd_->file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&Faet);
    global_dpd_->file2_mat_close(&Faet);
    global_dpd_->file2_close(&Faet);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_mat_init(&FAEt);
    global_dpd_->file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FAEt);
    global_dpd_->file2_mat_close(&FAEt);
    global_dpd_->file2_close(&FAEt);

    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 0, 0, "Fmit");
    global_dpd_->file2_mat_init(&Fmit);
    global_dpd_->file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&Fmit);
    global_dpd_->file2_mat_close(&Fmit);
    global_dpd_->file2_close(&Fmit);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_mat_init(&FMIt);
    global_dpd_->file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FMIt);
    global_dpd_->file2_mat_close(&FMIt);
    global_dpd_->file2_close(&FMIt);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* FME and Fme */
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->contract422(&D, &tia, &FME, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

    /* FAE and Fae */
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&fAB, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&fAB);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->file2_copy(&fab, PSIF_CC_OEI, "Fae");
    global_dpd_->file2_close(&fab);

    if((params.wfn == "CC2") || params.wfn == "EOM_CC2" ) {
      global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
      global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

      global_dpd_->file2_mat_init(&FAE);
      global_dpd_->file2_mat_rd(&FAE);
      global_dpd_->file2_mat_init(&Fae);
      global_dpd_->file2_mat_rd(&Fae);

      for(h=0; h < moinfo.nirreps; h++) {
	for(e=0; e < FAE.params->coltot[h]; e++)
	  FAE.matrix[h][e][e] = 0;
	for(e=0; e < Fae.params->coltot[h]; e++)
	  Fae.matrix[h][e][e] = 0;
      }

      global_dpd_->file2_mat_wrt(&FAE);
      global_dpd_->file2_mat_close(&FAE);
      global_dpd_->file2_mat_wrt(&Fae);
      global_dpd_->file2_mat_close(&Fae);

      global_dpd_->file2_close(&FAE);
      global_dpd_->file2_close(&Fae);
    }

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->dot13(&tia, &F, &Fae, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    global_dpd_->contract442(&taut, &D, &FAE, 3, 3, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    global_dpd_->contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    global_dpd_->buf4_close(&taut);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->contract222(&tia, &Fme, &Fae, 1, 1, -0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&Fme);

    /* Fae and FAE tilde intermediates */
    global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
    global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

    global_dpd_->file2_close(&FAE);
    global_dpd_->file2_close(&Fae);

    /* FMI and Fmi */
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_copy(&fIJ, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&fIJ);

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    global_dpd_->file2_copy(&fij, PSIF_CC_OEI, "Fmi");
    global_dpd_->file2_close(&fij);

    if( params.wfn == "CC2" || params.wfn == "EOM_CC2" ) {
      global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
      global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

      global_dpd_->file2_mat_init(&FMI);
      global_dpd_->file2_mat_rd(&FMI);
      global_dpd_->file2_mat_init(&Fmi);
      global_dpd_->file2_mat_rd(&Fmi);

      for(h=0; h < moinfo.nirreps; h++) {
	for(m=0; m < FMI.params->rowtot[h]; m++)
	  FMI.matrix[h][m][m] = 0;
	for(m=0; m < Fmi.params->rowtot[h]; m++)
	  Fmi.matrix[h][m][m] = 0;
      }

      global_dpd_->file2_mat_wrt(&FMI);
      global_dpd_->file2_mat_close(&FMI);
      global_dpd_->file2_mat_wrt(&Fmi);
      global_dpd_->file2_mat_close(&Fmi);

      global_dpd_->file2_close(&FMI);
      global_dpd_->file2_close(&Fmi);
    }

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&tia);
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

    global_dpd_->dot13(&tIA, &E_anti, &FMI, 1, 1, 1, 1);
    global_dpd_->dot24(&tia, &E, &FMI, 0, 0, 1, 1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&E_anti);

    global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

    global_dpd_->dot13(&tia, &E_anti, &Fmi, 1, 1, 1, 1);
    global_dpd_->dot13(&tIA, &E, &Fmi, 1, 1, 1, 1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&E_anti);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    global_dpd_->contract442(&D, &tautIJAB, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIJAB);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tautijab");
    global_dpd_->contract442(&D, &tautijab, &Fmi, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautijab);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    global_dpd_->contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tautiJaB");
    global_dpd_->contract442(&D, &tautIjAb, &Fmi, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tautIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&tIA);

    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->contract222(&Fme, &tia, &Fmi, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&Fme);
    global_dpd_->file2_close(&tia);

    /* FMI and Fmi tilde intermediate */
    global_dpd_->file2_copy(&FMI, PSIF_CC_OEI, "FMIt");
    global_dpd_->file2_copy(&Fmi, PSIF_CC_OEI, "Fmit");

    global_dpd_->file2_close(&FMI);
    global_dpd_->file2_close(&Fmi);

    /* remove diagonal elements from Ft's */
    global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 3, 3, "Faet");
    global_dpd_->file2_mat_init(&Faet);
    global_dpd_->file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&Faet);
    global_dpd_->file2_mat_close(&Faet);
    global_dpd_->file2_close(&Faet);

    global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
    global_dpd_->file2_mat_init(&FAEt);
    global_dpd_->file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FAEt);
    global_dpd_->file2_mat_close(&FAEt);
    global_dpd_->file2_close(&FAEt);

    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 2, 2, "Fmit");
    global_dpd_->file2_mat_init(&Fmit);
    global_dpd_->file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&Fmit);
    global_dpd_->file2_mat_close(&Fmit);
    global_dpd_->file2_close(&Fmit);

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_mat_init(&FMIt);
    global_dpd_->file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    global_dpd_->file2_mat_wrt(&FMIt);
    global_dpd_->file2_mat_close(&FMIt);
    global_dpd_->file2_close(&FMIt);
  } /** UHF **/

  return;

}

}} // namespace psi::cchbar
