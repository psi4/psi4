/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/*** This function removes incorrectly non-zero elements from ***
 *** a vector.  The non-zero elements are due to the          ***
 *** specification of open-shell orbitals as both occupied    ***
 *** and virtual orbitals.                                    ***/

void c_clean(dpdfile2 *CME, dpdfile2 *Cme,
 dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf) {

  int *occpi, *virtpi, *occ_off, *vir_off, *openpi, C_irr;
  int nirreps, *occ_sym, *vir_sym;
  int mn, ef, m, n, e, f, h, M, N, E, F;
  int msym, nsym, esym, fsym;

  C_irr = CME->my_irrep;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file2_mat_init(CME);
  dpd_file2_mat_rd(CME);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h^C_irr]-openpi[h^C_irr]); e<virtpi[h^C_irr]; e++)
        CME->matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(CME);

  dpd_file2_mat_init(Cme);
  dpd_file2_mat_rd(Cme);
  for(h=0; h < nirreps; h++) {
    for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
      for(e=0; e<virtpi[h^C_irr]; e++)
        Cme->matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(Cme);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(CMNEF, h);
    dpd_buf4_mat_irrep_rd(CMNEF, h);
    for(mn=0; mn < CMNEF->params->rowtot[h]; mn++) {
      for(ef=0; ef < CMNEF->params->coltot[h^C_irr]; ef++) {
          e = CMNEF->params->colorb[h^C_irr][ef][0];
          f = CMNEF->params->colorb[h^C_irr][ef][1];
          esym = CMNEF->params->rsym[e];
          fsym = CMNEF->params->ssym[f];
          E = e - vir_off[esym];
          F = f - vir_off[fsym];
          if ((E >= (virtpi[esym] - openpi[esym])) ||
              (F >= (virtpi[fsym] - openpi[fsym])) )
                   CMNEF->matrix[h][mn][ef] = 0.0;
      }
    }
    dpd_buf4_mat_irrep_wrt(CMNEF, h);
    dpd_buf4_mat_irrep_close(CMNEF, h);
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(Cmnef, h);
    dpd_buf4_mat_irrep_rd(Cmnef, h);
    for(mn=0; mn < Cmnef->params->rowtot[h]; mn++) {
      m = Cmnef->params->roworb[h][mn][0];
      n = Cmnef->params->roworb[h][mn][1];
      msym = Cmnef->params->psym[m];
      nsym = Cmnef->params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ef=0; ef < Cmnef->params->coltot[h^C_irr]; ef++) {
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (N >= (occpi[nsym] - openpi[nsym])) )
               Cmnef->matrix[h][mn][ef] = 0.0;
      }
    }
    dpd_buf4_mat_irrep_wrt(Cmnef, h);
    dpd_buf4_mat_irrep_close(Cmnef, h);
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(CMnEf, h);
    dpd_buf4_mat_irrep_rd(CMnEf, h);
    for(mn=0; mn < CMnEf->params->rowtot[h]; mn++) {
      n = CMnEf->params->roworb[h][mn][1];
      nsym = CMnEf->params->qsym[n];
      N = n - occ_off[nsym];
      for(ef=0; ef < CMnEf->params->coltot[h^C_irr]; ef++) {
        e = CMnEf->params->colorb[h^C_irr][ef][0];
        esym = CMnEf->params->rsym[e];
        E = e - vir_off[esym];
        if ((N >= (occpi[nsym] - openpi[nsym])) ||
            (E >= (virtpi[esym] - openpi[esym])) )
          CMnEf->matrix[h][mn][ef] = 0.0;
      }
    }
    dpd_buf4_mat_irrep_wrt(CMnEf, h);
    dpd_buf4_mat_irrep_close(CMnEf, h);
  }

  return;
}


void c_cleanSS(dpdfile2 *CME, dpdfile2 *Cme) {
  int *occpi, *virtpi, *occ_off, *vir_off, *openpi;
  int nirreps, *occ_sym, *vir_sym;
  int mn, ef, m, n, e, f;
  int h, M, N, E, F;
  int msym, nsym, esym, fsym, C_irr;

  C_irr = CME->my_irrep;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file2_mat_init(CME);
  dpd_file2_mat_rd(CME);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h^C_irr]-openpi[h^C_irr]); e<virtpi[h^C_irr]; e++)
        CME->matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(CME);

  dpd_file2_mat_init(Cme);
  dpd_file2_mat_rd(Cme);
  for(h=0; h < nirreps; h++) {
    for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
      for(e=0; e<virtpi[h^C_irr]; e++)
        Cme->matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(Cme);

  return;
}

}} // namespace psi::cceom
