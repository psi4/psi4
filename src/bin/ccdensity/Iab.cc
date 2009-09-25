/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* Iab(): Build the virtual-virtual block of the orbital Lagrangian
    ** using the expression given in lag.c.
    **
    ** Note that the code currently produces only the unique I_AB terms,
    ** but the actual contractions still need to be spin-adapted for
    ** greater efficiency.
    **
    ** TDC, 2/2008
    */

    void Iab(struct RHO_Params rho_params)
    {
      int a, b, c, A, B, C, Ga, Gb, Gc, Gac, Gbc;
      int *vir_off, *virtpi, nirreps, length, col;
      dpdfile2 F, D, I;
      dpdbuf4 G, Bints, Cints, Dints, Eints, Fints;
      vir_off = moinfo.vir_off;
      virtpi = moinfo.virtpi;
      nirreps = moinfo.nirreps;

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- sum_I fAI (DBI + DIB) + sum_C fAC (DBC + DCB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- sum_I fAI (DBI + DIB) + sum_C fAC (DBC + DCB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ab <-- sum_i fai (Dbi + Dib) + sum_c fac (Dbc + Dcb) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fab");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- sum_I fAI (DBI + DIB) + sum_C fAC (DBC + DCB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
	dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

	/* I'ab <-- sum_i fai (Dbi + Dib) + sum_c fac (Dbc + Dcb) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 0.0);
	dpd_file2_close(&D);
	dpd_file2_init(&D, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
	dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_init(&F, CC_OEI, 0, 3, 3, "fab");
	dpd_file2_init(&D, CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
	dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
	dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
	dpd_file2_close(&D);
	dpd_file2_close(&F);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- sum_JKI <JK||IA> G(JK,IB) + 2 sum_jKi <jK|iA> G(jK,iB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijka - Gjika", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "2 Gijka - Gjika");
	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&Eints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- sum_JKI <JK||IA> G(JK,IB) + 2 sum_jKi <jK|iA> G(jK,iB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_jki <jk||ia> G(jk,ib) + 2 sum_JkI <Jk|Ia> G(Jk,Ib) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- sum_JKI <JK||IA> G(JK,IB) + 2 sum_jKi <jK|iA> G(jK,iB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 20, 2, 20, 0, "GIJKA");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_jki <jk||ia> G(jk,ib) + 2 sum_JkI <Jk|Ia> G(Jk,Ib) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Eints, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 30, 12, 30, 0, "Gijka");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_buf4_init(&Eints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
	dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Eints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- sum_CDE <AC||DE> G(BC,DE) + 2 sum_cDe <Ac|De> G(Bc,De) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gabcd - Gabdc", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "2 Gabcd - Gabdc");
	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Bints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- sum_CDE <AC||DE> G(BC,DE) + 2 sum_cDe <Ac|De> G(Bc,De) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "GABCD");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	/*    dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);  replaced with 2(V**3) memory code*/
	dpd_file2_mat_init(&I);
	dpd_file2_mat_rd(&I);
	for(Gac=0; Gac < nirreps; Gac++) {
	  Gbc = Gac;
	  for(Ga=0; Ga < nirreps; Ga++) {
	    Gb = Ga;
	    Gc = Ga ^ Gac;
	    Bints.matrix[Gac] = dpd_block_matrix(virtpi[Gc], Bints.params->coltot[Gac]);
	    G.matrix[Gbc] = dpd_block_matrix(virtpi[Gc], G.params->coltot[Gbc]);
	    for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < virtpi[Gb]; b++) {
		B = vir_off[Gb] + b;
		dpd_buf4_mat_irrep_rd_block(&Bints, Gac, Bints.row_offset[Gac][A], virtpi[Gc]);
		dpd_buf4_mat_irrep_rd_block(&G, Gbc, G.row_offset[Gbc][B], virtpi[Gc]);
		length = virtpi[Gc] * Bints.params->coltot[Gac];
		if(length)
		  I.matrix[Ga][a][b] += 2.0 * C_DDOT(length, Bints.matrix[Gac][0], 1, G.matrix[Gbc][0], 1);
	      }
	    }
	    dpd_free_block(Bints.matrix[Gac], virtpi[Gc], Bints.params->coltot[Gac]);
	    dpd_free_block(G.matrix[Gbc], virtpi[Gc], G.params->coltot[Gbc]);
	  }
	}
	dpd_file2_mat_wrt(&I);
	dpd_file2_mat_close(&I);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_cde <ac||de> G(bc,de) + 2 sum_CdE <Ed|Ca> G(Ed,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "Gabcd");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	/*  dpd_contract442(&Bints, &G, &I, 3, 3, 2.0, 1.0); replaced with 2(V**3) memory code*/
	dpd_file2_mat_init(&I);
	dpd_file2_mat_rd(&I);
	for(Gac=0; Gac < nirreps; Gac++) {
	  Gbc = Gac;
	  for(Gc=0; Gc < nirreps; Gc++) {
	    Ga = Gc ^ Gac;
	    Gb = Ga;
	    Bints.matrix[Gac] = dpd_block_matrix(virtpi[Ga], Bints.params->coltot[Gac]);
	    G.matrix[Gbc] = dpd_block_matrix(virtpi[Gb], G.params->coltot[Gbc]);
	    for(c=0; c < virtpi[Gc]; c++) {
	      C = vir_off[Gc] + c;
	      dpd_buf4_mat_irrep_rd_block(&Bints, Gac, Bints.row_offset[Gac][C], virtpi[Ga]);
	      dpd_buf4_mat_irrep_rd_block(&G, Gbc, G.row_offset[Gbc][C], virtpi[Gb]);
	      for(a=0; a < virtpi[Ga]; a++) {
		for(b=0; b < virtpi[Gb]; b++) {
		  for (col=0; col< Bints.params->coltot[Gac]; ++col)
		    I.matrix[Ga][a][b] += 2.0 * Bints.matrix[Gac][a][col]*G.matrix[Gbc][b][col] ;
		}
	      }
	    }
	    dpd_free_block(Bints.matrix[Gac], virtpi[Ga], Bints.params->coltot[Gac]);
	    dpd_free_block(G.matrix[Gbc], virtpi[Gb], G.params->coltot[Gbc]);
	  }
	}
	dpd_file2_mat_wrt(&I);
	dpd_file2_mat_close(&I);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- sum_CDE <AC||DE> G(BC,DE) + 2 sum_cDe <Ac|De> G(Bc,De) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 7, 5, 5, 1, "B <AB|CD>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "GABCD");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_buf4_init(&Bints, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_cde <ac||de> G(bc,de) + 2 sum_CdE <Ed|Ca> G(Ed,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Bints, CC_BINTS, 0, 15, 17, 15, 15, 1, "B <ab|cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 15, 17, 17, 17, 0, "Gabcd");
	dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_buf4_init(&Bints, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 28, 28, 28, 28, 0, "GAbCd");
	dpd_contract442(&Bints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Bints);

	dpd_file2_close(&I);

      }

  
 
      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- sum_ICD <AI||CD> G(BI,CD) + 2 sum_iCd <Ai|Cd> G(Bi,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	dpd_buf4_init(&Fints, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&Fints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- sum_ICD <AI||CD> G(BI,CD) + 2 sum_iCd <Ai|Cd> G(Bi,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qprs, 11, 7, "F(CI,AB)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 7, 11, 7, 0, "F(CI,AB)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F(Ai,Cd)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F(Ai,Cd)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_icd <ai||cd> G(bi,cd) + 2 sum_IcD <aI|cD> G(bI,cD) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qprs, 11, 7, "F(ci,ab)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 7, 11, 7, 0, "F(ci,ab)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
	dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F(aI,cD)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F(aI,cD)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- sum_ICD <AI||CD> G(BI,CD) + 2 sum_iCd <Ai|Cd> G(Bi,Cd) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_icd <ai||cd> G(bi,cd) + 2 sum_IcD <aI|cD> G(bI,cD) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- 2 sum_CDI <DI||CA> G(DI,CB) + 2 sum_cDi <Di|Ac> G(Di,Bc)
	   + 2 sum_cdI <dI|cA> G(dI,cB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gciab - Gciba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "2 Gciab - Gciba");
	dpd_buf4_init(&Fints, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&Fints);
	dpd_buf4_close(&G);


	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- 2 sum_CDI <DI||CA> G(DI,CB) + 2 sum_cDi <Di|Ac> G(Di,Bc)
	   + 2 sum_cdI <dI|cA> G(dI,cB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qprs, 11, 5, "F (DI,CA)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F (DI,CA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
	dpd_contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_buf4_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F (Di,Ac)");
	dpd_buf4_close(&Fints);
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F (Di,Ac)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ab <-- 2 sum_cdi <di||ca> G(di,cb) + 2 sum_CdI <dI|aC> G(dI,bC)
	   + 2 sum_CDi <Di|Ca> G(Di,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	/* Both sorted F-blocks used here were generated above */
	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F (DI,CA)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
	dpd_contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_TMP0, 0, 11, 5, 11, 5, 0, "F (Di,Ac)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/ 

	/* I'AB <-- 2 sum_CDI <DI||CA> G(DI,CB) + 2 sum_cDi <Di|Ac> G(Di,Bc)
	   + 2 sum_cdI <dI|cA> G(dI,cB) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 21, 5, 21, 5, 1, "F <AI|BC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);

	/* I'ab <-- 2 sum_cdi <di||ca> G(di,cb) + 2 sum_CdI <dI|aC> G(dI,bC)
	   + 2 sum_CDi <Di|Ca> G(Di,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Fints, CC_FINTS, 0, 31, 15, 31, 15, 1, "F <ai|bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
	dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_buf4_init(&Fints, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
	dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Fints);

	dpd_file2_close(&I);
      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- 2 sum_IJC <JC||IA> G(JC,IB) + 2 sum_jCi <jC|iA> G(jC,iB)
	   - 2 sum_Jci <Jc|Ai> G(Jc,Bi) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, rpqs, 0, 5, "GIbjA (jI,bA)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (jI,bA)");
	dpd_contract442(&Dints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- 2 sum_IJC <JC||IA> G(JC,IB) + 2 sum_jCi <jC|iA> G(jC,iB)
	   - 2 sum_Jci <Jc|Ai> G(Jc,Bi) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, rpqs, 0, 5, "GIbjA (jI,bA)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (jI,bA)");
	dpd_contract442(&Dints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ab <-- 2 sum_jci <jc||ia> G(jc,ib) + 2 sum_JcI <Jc|Ia> G(Jc,Ib)
	   - 2 sum_jCI <Ij|Ca> GjCbI (IC,jb) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, rpqs, 0, 5, "GiBJa (Ji,Ba)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GiBJa (Ji,Ba)");
	dpd_contract442(&Dints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- 2 sum_IJC <JC||IA> G(JC,IB) + 2 sum_jCi <jC|iA> G(jC,iB)
	   - 2 sum_Jci <Jc|Ai> G(Jc,Bi) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
	dpd_buf4_sort(&G, CC_TMP0, rpqs, 23, 29, "GIbjA (jI,bA)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 23, 29, 23, 29, 0, "GIbjA (jI,bA)");
	dpd_contract442(&Dints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ab <-- 2 sum_jci <jc||ia> G(jc,ib) + 2 sum_JcI <Jc|Ia> G(Jc,Ib)
	   - 2 sum_jCI <Ij|Ca> GjCbI (IC,jb) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Cints, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Cints, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
	dpd_contract442(&Cints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Cints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
	dpd_buf4_sort(&G, CC_TMP0, rpqs, 22, 28, "GiBJa (Ji,Ba)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 22, 28, 22, 28, 0, "GiBJa (Ji,Ba)");
	dpd_contract442(&Dints, &G, &I, 3, 3, -2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }

  

      if(params.ref == 0) { /** RHF **/
	/* I'AB <-- sum_CIJ <IJ||CA> G(IJ,CB) + 2 sum_Ijc <Ij|Ac> G(Ij,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_buf4_scmcopy(&G, CC_GAMMA, "2 Gijab - Gijba", 2);
	dpd_buf4_sort_axpy(&G, CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&Dints);
	dpd_buf4_close(&G);

	dpd_file2_close(&I);
      }
      else if(params.ref == 1) { /** ROHF **/

	/* I'AB <-- sum_CIJ <IJ||CA> G(IJ,CB) + 2 sum_Ijc <Ij|Ac> G(Ij,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_cij <ij||ca> G(ij,cb) + 2 sum_IjC <Ij|Ca> G(Ij,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "Gijab");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);
      }
      else if(params.ref == 2) { /** UHF **/

	/* I'AB <-- sum_CIJ <IJ||CA> G(IJ,CB) + 2 sum_Ijc <Ij|Ac> G(Ij,Bc) */
	dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

	/* I'ab <-- sum_cij <ij||ca> G(ij,cb) + 2 sum_IjC <Ij|Ca> G(Ij,Cb) */
	dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");

	dpd_buf4_init(&Dints, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 15, 12, 17, 0, "Gijab");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_buf4_init(&Dints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_contract442(&Dints, &G, &I, 3, 3, 2.0, 1.0);
	dpd_buf4_close(&G);
	dpd_buf4_close(&Dints);

	dpd_file2_close(&I);

      }

    }

  }} // namespace psi::ccdensity
