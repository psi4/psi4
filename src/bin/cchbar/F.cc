/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/

/*! \defgroup CCHBAR cchbar: Compute the similarity-transformed Hamiltonian */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
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
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);

    /** FAE **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");

    /*
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F 2<ia|bc> - <ia|cb>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&F);
    */

    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&tIA);
    dpd_file2_mat_rd(&tIA);
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    for(Gma=0; Gma < moinfo.nirreps; Gma++) {
      dpd_buf4_mat_irrep_row_init(&F, Gma);
      X = init_array(F.params->coltot[Gma]);

      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	dpd_buf4_mat_irrep_row_rd(&F, Gma, ma);
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
      dpd_buf4_mat_irrep_row_close(&F, Gma);
    }
    dpd_buf4_close(&F);
    dpd_file2_mat_close(&tIA);
    dpd_file2_close(&tIA);
    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&FAE);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_contract442(&tautIjAb, &D, &FAE, 3, 3, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_close(&FAE);

    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_dot13(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&E);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_close(&FMI);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_mat_init(&FAEt);
    dpd_file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FAEt);
    dpd_file2_mat_close(&FAEt);
    dpd_file2_close(&FAEt);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_mat_init(&FMIt);
    dpd_file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FMIt);
    dpd_file2_mat_close(&FMIt);
    dpd_file2_close(&FMIt);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* FME and Fme */
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fme");
    dpd_file2_close(&fia);
  
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  
    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);

    /* FAE and Fae */

    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
    dpd_file2_copy(&fab, CC_OEI, "Fae");
    dpd_file2_close(&fab);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    dpd_dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&F_anti);
    dpd_buf4_close(&F);

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    dpd_contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    dpd_contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&tautijab);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    dpd_contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&tia, &Fme, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);

    /* Form Fae and FAE tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);


    /* FMI and Fmi */

    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);
  
    dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
    dpd_file2_copy(&fij, CC_OEI, "Fmi");
    dpd_file2_close(&fij);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
  
    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);
  
    dpd_buf4_init(&E_anti, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
    dpd_dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0);

    dpd_dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0);
    dpd_dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&E_anti);
    dpd_buf4_close(&E);

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautijab");

    dpd_contract442(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_contract442(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1);

    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&tautijab);
    dpd_buf4_close(&D_anti);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_contract442(&D, &tautIjAb, &Fmi, 1, 1, 1, 1);

    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&Fme, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&Fme);
    dpd_file2_close(&tia);

    /* FMI and Fmi tilde intermediates */
    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_copy(&Fmi, CC_OEI, "Fmit");

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);


    /* remove diagonal elements from Ft's */
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");
    dpd_file2_mat_init(&Faet);
    dpd_file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Faet);
    dpd_file2_mat_close(&Faet);
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_mat_init(&FAEt);
    dpd_file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FAEt);
    dpd_file2_mat_close(&FAEt);
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 0, 0, "Fmit");
    dpd_file2_mat_init(&Fmit);
    dpd_file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Fmit);
    dpd_file2_mat_close(&Fmit);
    dpd_file2_close(&Fmit);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_mat_init(&FMIt);
    dpd_file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FMIt);
    dpd_file2_mat_close(&FMIt);
    dpd_file2_close(&FMIt);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* FME and Fme */
    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_copy(&fIA, CC_OEI, "FME");
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_copy(&fia, CC_OEI, "Fme");
    dpd_file2_close(&fia);
  
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_contract422(&D, &tia, &FME, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);

    /* FAE and Fae */
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
    dpd_file2_copy(&fab, CC_OEI, "Fae");
    dpd_file2_close(&fab);

    if((!strcmp(params.wfn,"CC2")) || (!strcmp(params.wfn,"EOM_CC2"))) {
      dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
      dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");

      dpd_file2_mat_init(&FAE);
      dpd_file2_mat_rd(&FAE);
      dpd_file2_mat_init(&Fae);
      dpd_file2_mat_rd(&Fae);
  
      for(h=0; h < moinfo.nirreps; h++) {
	for(e=0; e < FAE.params->coltot[h]; e++)
	  FAE.matrix[h][e][e] = 0;
	for(e=0; e < Fae.params->coltot[h]; e++)
	  Fae.matrix[h][e][e] = 0;
      }

      dpd_file2_mat_wrt(&FAE);
      dpd_file2_mat_close(&FAE);
      dpd_file2_mat_wrt(&Fae);
      dpd_file2_mat_close(&Fae);

      dpd_file2_close(&FAE);
      dpd_file2_close(&Fae);
    }

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_dot13(&tia, &F, &Fae, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_contract442(&taut, &D, &FAE, 3, 3, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    dpd_contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_contract222(&tia, &Fme, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);

    /* Fae and FAE tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae); 

    /* FMI and Fmi */
    dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_copy(&fIJ, CC_OEI, "FMI");
    dpd_file2_close(&fIJ);
  
    dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
    dpd_file2_copy(&fij, CC_OEI, "Fmi");
    dpd_file2_close(&fij);

    if((!strcmp(params.wfn,"CC2")) || (!strcmp(params.wfn,"EOM_CC2"))) {
      dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
      dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");

      dpd_file2_mat_init(&FMI);
      dpd_file2_mat_rd(&FMI);
      dpd_file2_mat_init(&Fmi);
      dpd_file2_mat_rd(&Fmi);

      for(h=0; h < moinfo.nirreps; h++) {
	for(m=0; m < FMI.params->rowtot[h]; m++)
	  FMI.matrix[h][m][m] = 0;
	for(m=0; m < Fmi.params->rowtot[h]; m++)
	  Fmi.matrix[h][m][m] = 0;
      }

      dpd_file2_mat_wrt(&FMI);
      dpd_file2_mat_close(&FMI);
      dpd_file2_mat_wrt(&Fmi);
      dpd_file2_mat_close(&Fmi);

      dpd_file2_close(&FMI);
      dpd_file2_close(&Fmi);
    }

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);
  
    dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&fia, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&E_anti, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");

    dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1, 1);
    dpd_dot24(&tia, &E, &FMI, 0, 0, 1, 1);

    dpd_buf4_close(&E);
    dpd_buf4_close(&E_anti);

    dpd_buf4_init(&E_anti, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_init(&E, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");

    dpd_dot13(&tia, &E_anti, &Fmi, 1, 1, 1, 1);
    dpd_dot13(&tIA, &E, &Fmi, 1, 1, 1, 1);

    dpd_buf4_close(&E);
    dpd_buf4_close(&E_anti);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tautIJAB");
    dpd_contract442(&D, &tautIJAB, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tautijab");
    dpd_contract442(&D, &tautijab, &Fmi, 0, 0, 1, 1);
    dpd_buf4_close(&tautijab);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tautiJaB");
    dpd_contract442(&D, &tautIjAb, &Fmi, 0, 0, 1, 1);
    dpd_buf4_close(&tautIjAb);
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&FME, &tIA, &FMI, 0, 0, 0.5, 1);
    dpd_file2_close(&FME);
    dpd_file2_close(&tIA);

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_contract222(&Fme, &tia, &Fmi, 0, 0, 0.5, 1);
    dpd_file2_close(&Fme);
    dpd_file2_close(&tia);

    /* FMI and Fmi tilde intermediate */
    dpd_file2_copy(&FMI, CC_OEI, "FMIt");
    dpd_file2_copy(&Fmi, CC_OEI, "Fmit");

    dpd_file2_close(&FMI);
    dpd_file2_close(&Fmi);

    /* remove diagonal elements from Ft's */
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");
    dpd_file2_mat_init(&Faet);
    dpd_file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Faet);
    dpd_file2_mat_close(&Faet);
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_mat_init(&FAEt);
    dpd_file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FAEt);
    dpd_file2_mat_close(&FAEt);
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 2, 2, "Fmit");
    dpd_file2_mat_init(&Fmit);
    dpd_file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Fmit);
    dpd_file2_mat_close(&Fmit);
    dpd_file2_close(&Fmit);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_mat_init(&FMIt);
    dpd_file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FMIt);
    dpd_file2_mat_close(&FMIt);
    dpd_file2_close(&FMIt);
  } /** UHF **/

  return;

}

}} // namespace psi::cchbar
