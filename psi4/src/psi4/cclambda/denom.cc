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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void denom_rhf(struct L_Params);
void denom_rohf(struct L_Params);
void denom_uhf(struct L_Params);

void denom(struct L_Params L_params) {
  if(params.ref == 0) denom_rhf(L_params);
  else if(params.ref == 1) denom_rohf(L_params);
  else if(params.ref == 2) denom_uhf(L_params);
}

void denom_rhf(struct L_Params L_params)
{
  dpdfile2 FAE, FMI;
  dpdfile2 dIA;
  dpdfile4 dIjAb;
  dpdbuf4 d, bdIJAB, bdijab, bdIjAb;
  double tval;
  int nirreps,L_irr;
  int h, i, j, a, b, ij, ab;
  int I, J, A, B;
  int isym, jsym, asym, bsym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *openpi;
  double Fii, Fjj, Faa, Fbb;

  L_irr = L_params.irrep;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
  global_dpd_->file2_mat_init(&FMI);
  global_dpd_->file2_mat_rd(&FMI);

  global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
  global_dpd_->file2_mat_init(&FAE);
  global_dpd_->file2_mat_rd(&FAE);

  /* Alpha one-electron denominator */
  global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) { /* irreps of dIA and Fii */
    for(i=0; i < occpi[h]; i++) {
      Fii = FMI.matrix[h][i][i];
      for(a=0; a < virtpi[h^L_irr]; a++) {
        Faa = FAE.matrix[h^L_irr][a][a];
        dIA.matrix[h][i][a] = 1.0/(Fii - Faa + L_params.cceom_energy);
      }
    }
  }
  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  /* Alpha-beta two-electron denominator */
  global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, "dIjAb");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dIjAb, h);
    /* Loop over the rows */
    for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
          i = dIjAb.params->roworb[h][ij][0];
          j = dIjAb.params->roworb[h][ij][1];
          isym = dIjAb.params->psym[i];
          jsym = dIjAb.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];
          Fii = FMI.matrix[isym][I][I];
          Fjj = FMI.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dIjAb.params->coltot[h^L_irr]; ab++) {
              a = dIjAb.params->colorb[h^L_irr][ab][0];
              b = dIjAb.params->colorb[h^L_irr][ab][1];
              asym = dIjAb.params->rsym[a];
              bsym = dIjAb.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = FAE.matrix[asym][A][A];
              Fbb = FAE.matrix[bsym][B][B];

              dIjAb.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb + L_params.cceom_energy);
            }
        }
    global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
    global_dpd_->file4_mat_irrep_close(&dIjAb, h);
  }
  global_dpd_->file4_close(&dIjAb);

  global_dpd_->file2_mat_close(&FMI);
  global_dpd_->file2_mat_close(&FAE);
  global_dpd_->file2_close(&FMI);
  global_dpd_->file2_close(&FAE);

  return;
}

void denom_uhf(struct L_Params L_params)
{
  int nirreps, h, i, j, a, b, ij, ab, I, J, A, B, isym, jsym, asym, bsym, m, e;
  int *aoccpi, *boccpi, *avirtpi, *bvirtpi;
  int *aocc_off, *bocc_off, *avir_off, *bvir_off, L_irr;
  dpdfile2 LFMIt, LFmit, LFaet, LFAEt;
  dpdfile2 FMI, Fmi, FAE, Fae;
  dpdfile2 dIA, dia;
  dpdfile4 dIJAB, dijab, dIjAb;
  double Fii, Fjj, Faa, Fbb;

  L_irr = L_params.irrep;
  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  boccpi = moinfo.boccpi;
  avirtpi = moinfo.avirtpi;
  bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off;
  bocc_off = moinfo.bocc_off;
  avir_off = moinfo.avir_off;
  bvir_off = moinfo.bvir_off;

  if((params.wfn == "CC2") || (params.wfn == "EOM_CC2")) {

    global_dpd_->file2_init(&LFMIt, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&LFMIt);
    global_dpd_->file2_mat_rd(&LFMIt);

    global_dpd_->file2_init(&LFmit, PSIF_CC_OEI, 0, 2, 2, "fij");
    global_dpd_->file2_mat_init(&LFmit);
    global_dpd_->file2_mat_rd(&LFmit);

    global_dpd_->file2_init(&LFaet, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->file2_mat_init(&LFaet);
    global_dpd_->file2_mat_rd(&LFaet);

    global_dpd_->file2_init(&LFAEt, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&LFAEt);
    global_dpd_->file2_mat_rd(&LFAEt);

  }
  else {
    global_dpd_->file2_init(&LFMIt, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_mat_init(&LFMIt);
    global_dpd_->file2_mat_rd(&LFMIt);

    global_dpd_->file2_init(&LFmit, PSIF_CC_OEI, 0, 2, 2, "Fmi");
    global_dpd_->file2_mat_init(&LFmit);
    global_dpd_->file2_mat_rd(&LFmit);

    global_dpd_->file2_init(&LFaet, PSIF_CC_OEI, 0, 3, 3, "Fae");
    global_dpd_->file2_mat_init(&LFaet);
    global_dpd_->file2_mat_rd(&LFaet);

    global_dpd_->file2_init(&LFAEt, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_mat_init(&LFAEt);
    global_dpd_->file2_mat_rd(&LFAEt);
  }

  global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      Fii = LFMIt.matrix[h][i][i];
      for(a=0; a < avirtpi[h^L_irr]; a++) {
        Faa = LFAEt.matrix[h^L_irr][a][a];
        dIA.matrix[h][i][a] = 1.0/(Fii - Faa + L_params.cceom_energy);
      }
    }
  }
  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  global_dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 2, 3, "dia");
  global_dpd_->file2_mat_init(&dia);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      Fii = LFmit.matrix[h][i][i];
      for(a=0; a < bvirtpi[h^L_irr]; a++) {
        Faa = LFaet.matrix[h^L_irr][a][a];
        dia.matrix[h][i][a] = 1.0/(Fii - Faa + L_params.cceom_energy);
      }
    }
  }
  global_dpd_->file2_mat_wrt(&dia);
  global_dpd_->file2_mat_close(&dia);
  global_dpd_->file2_close(&dia);

  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 1, 6, "dIJAB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dIJAB, h);
    for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
      i = dIJAB.params->roworb[h][ij][0];
      j = dIJAB.params->roworb[h][ij][1];
      isym = dIJAB.params->psym[i];
      jsym = dIJAB.params->qsym[j];
      I = i - aocc_off[isym];
      J = j - aocc_off[jsym];
      Fii = LFMIt.matrix[isym][I][I];
      Fjj = LFMIt.matrix[jsym][J][J];

      for(ab=0; ab < dIJAB.params->coltot[h^L_irr]; ab++) {
	a = dIJAB.params->colorb[h^L_irr][ab][0];
	b = dIJAB.params->colorb[h^L_irr][ab][1];
	asym = dIJAB.params->rsym[a];
	bsym = dIJAB.params->ssym[b];
	A = a - avir_off[asym];
	B = b - avir_off[bsym];
	Faa = LFAEt.matrix[asym][A][A];
	Fbb = LFAEt.matrix[bsym][B][B];

	dIJAB.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
				       + L_params.cceom_energy);
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
    global_dpd_->file4_mat_irrep_close(&dIJAB, h);
  }
  global_dpd_->file4_close(&dIJAB);

  global_dpd_->file4_init(&dijab, PSIF_CC_DENOM, L_irr, 11, 16, "dijab");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dijab, h);
    for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
      i = dijab.params->roworb[h][ij][0];
      j = dijab.params->roworb[h][ij][1];
      isym = dijab.params->psym[i];
      jsym = dijab.params->qsym[j];
      I = i - bocc_off[isym];
      J = j - bocc_off[jsym];
      Fii = LFmit.matrix[isym][I][I];
      Fjj = LFmit.matrix[jsym][J][J];

      for(ab=0; ab < dijab.params->coltot[h^L_irr]; ab++) {
	a = dijab.params->colorb[h^L_irr][ab][0];
	b = dijab.params->colorb[h^L_irr][ab][1];
	asym = dijab.params->rsym[a];
	bsym = dijab.params->ssym[b];
	A = a - bvir_off[asym];
	B = b - bvir_off[bsym];
	Faa = LFaet.matrix[asym][A][A];
	Fbb = LFaet.matrix[bsym][B][B];

	dijab.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
				       + L_params.cceom_energy);
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&dijab, h);
    global_dpd_->file4_mat_irrep_close(&dijab, h);
  }
  global_dpd_->file4_close(&dijab);

  global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 22, 28, "dIjAb");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dIjAb, h);
    for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
      i = dIjAb.params->roworb[h][ij][0];
      j = dIjAb.params->roworb[h][ij][1];
      isym = dIjAb.params->psym[i];
      jsym = dIjAb.params->qsym[j];
      I = i - aocc_off[isym];
      J = j - bocc_off[jsym];
      Fii = LFMIt.matrix[isym][I][I];
      Fjj = LFmit.matrix[jsym][J][J];

      for(ab=0; ab < dIjAb.params->coltot[h^L_irr]; ab++) {
	a = dIjAb.params->colorb[h^L_irr][ab][0];
	b = dIjAb.params->colorb[h^L_irr][ab][1];
	asym = dIjAb.params->rsym[a];
	bsym = dIjAb.params->ssym[b];
	A = a - avir_off[asym];
	B = b - bvir_off[bsym];
	Faa = LFAEt.matrix[asym][A][A];
	Fbb = LFaet.matrix[bsym][B][B];

	dIjAb.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
				       + L_params.cceom_energy);
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
    global_dpd_->file4_mat_irrep_close(&dIjAb, h);
  }
  global_dpd_->file4_close(&dIjAb);

  global_dpd_->file2_mat_close(&LFMIt);
  global_dpd_->file2_mat_close(&LFmit);
  global_dpd_->file2_mat_close(&LFAEt);
  global_dpd_->file2_mat_close(&LFaet);
  global_dpd_->file2_close(&LFMIt);
  global_dpd_->file2_close(&LFmit);
  global_dpd_->file2_close(&LFAEt);
  global_dpd_->file2_close(&LFaet);

  /*     if((!strcmp(params.wfn,"CC2")) || (!strcmp(params.wfn,"EOM_CC2"))) { */
  /*       dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI"); */
  /*       dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi"); */

  /*       dpd_file2_mat_init(&FMI); */
  /*       dpd_file2_mat_rd(&FMI); */
  /*       dpd_file2_mat_init(&Fmi); */
  /*       dpd_file2_mat_rd(&Fmi); */

  /*       for(h=0; h < moinfo.nirreps; h++) { */
  /* 	for(m=0; m < FMI.params->rowtot[h]; m++)  */
  /* 	  FMI.matrix[h][m][m] = 0; */
  /* 	for(m=0; m < Fmi.params->rowtot[h]; m++)  */
  /* 	  Fmi.matrix[h][m][m] = 0; */
  /*       } */

  /*       dpd_file2_mat_wrt(&FMI); */
  /*       dpd_file2_mat_close(&FMI); */
  /*       dpd_file2_mat_wrt(&Fmi); */
  /*       dpd_file2_mat_close(&Fmi); */

  /*       dpd_file2_close(&FMI); */
  /*       dpd_file2_close(&Fmi); */

  /*       dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE"); */
  /*       dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae"); */

  /*       dpd_file2_mat_init(&FAE); */
  /*       dpd_file2_mat_rd(&FAE); */
  /*       dpd_file2_mat_init(&Fae); */
  /*       dpd_file2_mat_rd(&Fae); */

  /*       for(h=0; h < moinfo.nirreps; h++) { */
  /* 	for(e=0; e < FAE.params->coltot[h]; e++)  */
  /* 	  FAE.matrix[h][e][e] = 0; */
  /* 	for(e=0; e < Fae.params->coltot[h]; e++) */
  /* 	  Fae.matrix[h][e][e] = 0; */
  /*       } */

  /*       dpd_file2_mat_wrt(&FAE); */
  /*       dpd_file2_mat_close(&FAE); */
  /*       dpd_file2_mat_wrt(&Fae); */
  /*       dpd_file2_mat_close(&Fae); */

  /*       dpd_file2_close(&FAE); */
  /*       dpd_file2_close(&Fae); */
  /*     } */

  return;
}

void denom_rohf(struct L_Params L_params)
{
  dpdfile2 LFAEt, LFaet, LFMIt, LFmit;
  dpdfile2 dIA, dia;
  dpdfile4 dIJAB, dijab, dIjAb;
  dpdbuf4 d, bdIJAB, bdijab, bdIjAb;
  double tval;
  int nirreps,L_irr;
  int h, i, j, a, b, ij, ab;
  int I, J, A, B;
  int isym, jsym, asym, bsym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *openpi;
  double Fii, Fjj, Faa, Fbb;

  L_irr = L_params.irrep;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  global_dpd_->file2_init(&LFMIt, PSIF_CC_OEI, 0, 0, 0, "FMI");
  global_dpd_->file2_mat_init(&LFMIt);
  global_dpd_->file2_mat_rd(&LFMIt);

  global_dpd_->file2_init(&LFmit, PSIF_CC_OEI, 0, 0, 0, "Fmi");
  global_dpd_->file2_mat_init(&LFmit);
  global_dpd_->file2_mat_rd(&LFmit);

  global_dpd_->file2_init(&LFaet, PSIF_CC_OEI, 0, 1, 1, "Fae");
  global_dpd_->file2_mat_init(&LFaet);
  global_dpd_->file2_mat_rd(&LFaet);

  global_dpd_->file2_init(&LFAEt, PSIF_CC_OEI, 0, 1, 1, "FAE");
  global_dpd_->file2_mat_init(&LFAEt);
  global_dpd_->file2_mat_rd(&LFAEt);

  /* Alpha one-electron denominator */
  global_dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) { /* irreps of dIA and Fii */
    for(i=0; i < occpi[h]; i++) {
      Fii = LFMIt.matrix[h][i][i];
      for(a=0; a < (virtpi[h^L_irr] - openpi[h^L_irr]); a++) {
        Faa = LFAEt.matrix[h^L_irr][a][a];
        dIA.matrix[h][i][a] = 1.0/(Fii - Faa + L_params.cceom_energy);
      }
    }
  }
  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  /* Beta one-electron denominator */
  global_dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 0, 1, "dia");
  global_dpd_->file2_mat_init(&dia);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h] - openpi[h]); i++) {
      Fii = LFmit.matrix[h][i][i];
      for(a=0; a < virtpi[h^L_irr]; a++) {
        Faa = LFaet.matrix[h^L_irr][a][a];
        dia.matrix[h][i][a] = 1.0/(Fii - Faa + L_params.cceom_energy);
      }
    }
  }
  global_dpd_->file2_mat_wrt(&dia);
  global_dpd_->file2_mat_close(&dia);
  global_dpd_->file2_close(&dia);

  /* Alpha-alpha two-electron denominator */
  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 1, 6, "dIJAB");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dIJAB, h);
      /* Loop over the rows */
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
          i = dIJAB.params->roworb[h][ij][0];
          j = dIJAB.params->roworb[h][ij][1];
          isym = dIJAB.params->psym[i];
          jsym = dIJAB.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];

          Fii = LFMIt.matrix[isym][I][I];
          Fjj = LFMIt.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dIJAB.params->coltot[h^L_irr]; ab++) {
              a = dIJAB.params->colorb[h^L_irr][ab][0];
              b = dIJAB.params->colorb[h^L_irr][ab][1];
              asym = dIJAB.params->rsym[a];
              bsym = dIJAB.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFAEt.matrix[asym][A][A];
              Fbb = LFAEt.matrix[bsym][B][B];

              dIJAB.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
                 (B >= (virtpi[bsym] - openpi[bsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + L_params.cceom_energy));
          }
      }
    global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
    global_dpd_->file4_mat_irrep_close(&dIJAB, h);
  }
  global_dpd_->file4_close(&dIJAB);

  /* Beta-beta two-electron denominator */
  global_dpd_->file4_init(&dijab, PSIF_CC_DENOM, L_irr, 1, 6, "dijab");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dijab, h);
    /* Loop over the rows */
    for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
          i = dijab.params->roworb[h][ij][0];
          j = dijab.params->roworb[h][ij][1];
          isym = dijab.params->psym[i];
          jsym = dijab.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];

          Fii = LFmit.matrix[isym][I][I];
          Fjj = LFmit.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dijab.params->coltot[h^L_irr]; ab++) {
              a = dijab.params->colorb[h^L_irr][ab][0];
              b = dijab.params->colorb[h^L_irr][ab][1];
              asym = dijab.params->rsym[a];
              bsym = dijab.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFaet.matrix[asym][A][A];
              Fbb = LFaet.matrix[bsym][B][B];

              dijab.matrix[h][ij][ab] =
                ((I >= (occpi[isym] - openpi[isym])) ||
                 (J >= (occpi[jsym] - openpi[jsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + L_params.cceom_energy));
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&dijab, h);
    global_dpd_->file4_mat_irrep_close(&dijab, h);
  }
  global_dpd_->file4_close(&dijab);


  /* Alpha-beta two-electron denominator */
  global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, "dIjAb");

  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&dIjAb, h);
    /* Loop over the rows */
    for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
          i = dIjAb.params->roworb[h][ij][0];
          j = dIjAb.params->roworb[h][ij][1];
          isym = dIjAb.params->psym[i];
          jsym = dIjAb.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];
          Fii = LFMIt.matrix[isym][I][I];
          Fjj = LFmit.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dIjAb.params->coltot[h^L_irr]; ab++) {
              a = dIjAb.params->colorb[h^L_irr][ab][0];
              b = dIjAb.params->colorb[h^L_irr][ab][1];
              asym = dIjAb.params->rsym[a];
              bsym = dIjAb.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFAEt.matrix[asym][A][A];
              Fbb = LFaet.matrix[bsym][B][B];

              dIjAb.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
                 (J >= (occpi[jsym] - openpi[jsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + L_params.cceom_energy));
            }
        }
    global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
    global_dpd_->file4_mat_irrep_close(&dIjAb, h);
  }
  global_dpd_->file4_close(&dIjAb);

  global_dpd_->file2_mat_close(&LFMIt);
  global_dpd_->file2_mat_close(&LFmit);
  global_dpd_->file2_mat_close(&LFAEt);
  global_dpd_->file2_mat_close(&LFaet);
  global_dpd_->file2_close(&LFMIt);
  global_dpd_->file2_close(&LFmit);
  global_dpd_->file2_close(&LFAEt);
  global_dpd_->file2_close(&LFaet);

  return;
}




}} // namespace psi::cclambda
