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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void cc2_sigmaSS(int i, int C_irr);

void cc2_sigma(int i, int C_irr)
{
  dpdfile2 SIA;
  dpdbuf4 SIjAb;
  dpdfile2 CME;
  dpdbuf4 CMnEf;
  dpdfile2 FAE;
  dpdfile2 FMI;
  dpdfile2 FME;
  dpdbuf4 W;
  dpdbuf4 C;
  dpdfile2 S;
  dpdbuf4 WMbEj;
  dpdbuf4 WAmEf;
  dpdbuf4 WMnIe;
  dpdbuf4 WAbEi;
  dpdbuf4 WMbIj;
  dpdbuf4 Z;
  dpdbuf4 Z2;
  char lbl[32];
  char CME_lbl[32];
  char CMnEf_lbl[32];
  char SIjAb_lbl[32];
  int Gej, Gab, Gij, Gj, Gi, Ge, nrows, length, E, e, I;
  int Gam, Gef, Gim, Ga, Gm, ncols, A, a, am;

  if (params.eom_ref == 0) { /* RHF */

    cc2_sigmaSS(i,C_irr);

    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    global_dpd_->dot24(&FME,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&SIA);

    /*
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_contract442(&CMnEf, &WAmEf, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&SIA);
    */

    global_dpd_->buf4_init(&C, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&S, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    global_dpd_->file2_mat_init(&S);
    global_dpd_->file2_mat_rd(&S);
    for(Gam=0; Gam < moinfo.nirreps; Gam++) {
      Gef = Gam ^ H_IRR;
      Gim = Gef ^ C_irr;

      global_dpd_->buf4_mat_irrep_init(&C, Gim);
      global_dpd_->buf4_mat_irrep_rd(&C, Gim);
      global_dpd_->buf4_mat_irrep_shift13(&C, Gim);

      for(Gi=0; Gi < moinfo.nirreps; Gi++) {
        Ga = Gi ^ C_irr;
        Gm = Ga ^ Gam;

        W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gef]);

        nrows = moinfo.occpi[Gi];
        ncols = moinfo.occpi[Gm] * W.params->coltot[Gef];

        for(A=0; A < moinfo.virtpi[Ga]; A++) {
          a = moinfo.vir_off[Ga] + A;
          am = W.row_offset[Gam][a];

          global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

          if(nrows && ncols && moinfo.virtpi[Ga])
            C_DGEMV('n',nrows,ncols,1,C.shift.matrix[Gim][Gi][0],ncols,W.matrix[Gam][0], 1,
                    1, &(S.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
        }

        global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gef]);
      }

      global_dpd_->buf4_mat_irrep_close(&C, Gim);
    }
    global_dpd_->file2_mat_wrt(&S);
    global_dpd_->file2_mat_close(&S);
    global_dpd_->file2_close(&S);
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_close(&W);

    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
    global_dpd_->contract442(&WMnIe, &CMnEf, &SIA, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->file2_close(&SIA);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(Ij,Ab)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC2_HET1, H_IRR, 10, 0, 10, 0, 0, "CC2 WMbIj");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMbIj, &Z, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "WmaijDS Z(jI,bA)");
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb,  -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(jI,bA)");
    global_dpd_->buf4_axpy(&Z, &SIjAb,  -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(Ij,Ab)");
    global_dpd_->buf4_scm(&Z,0);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, H_IRR, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    /*dpd_contract244(&CME, &WAbEi, &Z, 1, 0, 0, 1.0, 0.0);*/
    global_dpd_->file2_mat_init(&CME);
    global_dpd_->file2_mat_rd(&CME);
    for(Gej=0; Gej < moinfo.nirreps; Gej++) {
      Gab = Gej ^ H_IRR;
      Gij = Gab ^ C_irr;

      global_dpd_->buf4_mat_irrep_init(&Z, Gij);
      global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);

      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
        Gj = Ge ^ Gej;
        Gi = Gj ^ Gij;

        nrows = moinfo.occpi[Gj];
        length = nrows * W.params->coltot[Gab];
        global_dpd_->buf4_mat_irrep_init_block(&W, Gej, nrows);

        for(E=0; E < moinfo.virtpi[Ge]; E++) {
          e = moinfo.vir_off[Ge] + E;
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][e], nrows);

          for(I=0; I < moinfo.occpi[Gi]; I++) {
            if(length)
              C_DAXPY(length, CME.matrix[Gi][I][E], W.matrix[Gej][0], 1,
                      Z.shift.matrix[Gij][Gi][I], 1);
          }
        }

        global_dpd_->buf4_mat_irrep_close_block(&W, Gej, nrows);
      }

      global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
      global_dpd_->buf4_mat_irrep_close(&Z, Gij);

    }
    global_dpd_->file2_mat_close(&CME);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&W);

    /*
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WabejDS Z(jI,bA)");
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
    */

    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_SIjAb, qpsr, 0, 5, SIjAb_lbl, 1);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1.0);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&Z);

    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(Ij,Ab)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "fAB");
    global_dpd_->contract424(&CMnEf, &FAE, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&FAE);
    global_dpd_->buf4_close(&CMnEf);

    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "FDD_Fbe Z(jI,bA)");
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(jI,bA)");

    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1.0);
    global_dpd_->buf4_axpy(&Z2, &SIjAb, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&SIjAb);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(Ij,Ab)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "fIJ");
    global_dpd_->contract244(&FMI, &CMnEf, &Z, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&FMI);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "FDD_Fmj Z(jI,bA)");
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(jI,bA)");
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);
  }
  else if (params.eom_ref == 1) { /* ROHF */
    printf("ROHF EOM_CC2 is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
  else { /* UHF */
    printf("UHF EOM_CC2 is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
}

void cc2_sigmaSS(int i, int C_irr)
{
  dpdfile2 SIA;
  dpdfile2 CME;
  dpdfile2 FAE;
  dpdfile2 FMI;
  dpdbuf4 WMbEj;
  dpdfile2 Xme;
  dpdbuf4 T2;
  dpdbuf4 D;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */

    /* sigmaSS */

    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
    global_dpd_->contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
    global_dpd_->contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&FMI);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC2_HET1, H_IRR, 10, 10, 10, 10, 0, "CC2 2 W(jb,ME) + W(Jb,Me)");
    global_dpd_->contract422(&WMbEj, &CME, &SIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->file2_init(&Xme, PSIF_CC_OEI, C_irr, 0, 1, "XME");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &CME, &Xme, 0, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    global_dpd_->contract422(&T2, &Xme, &SIA, 0, 0, 1, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&Xme);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&SIA);
  }
  else if (params.eom_ref == 1) { /* ROHF */
    printf("ROHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
  else { /* UHF */
    printf("UHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
}

}} // namespace psi::cceom
