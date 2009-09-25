/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* DUMP_ROHF(): Mulliken-order the ROHF-CCSD two-electron density and
** dump it to a file for subsequent backtransformation.  Basically all
** we have to do is swap indices two and three, e.g.
**
** G'(pr,qs) = G(pq,rs)
**
** In order for the Mulliken-ordered density to be valid for the
** backtransformation algorithm used in TRANSQT, the final density
** must have eight-fold permutational symmetry like the original
** integrals.  Unfortunately, there are a couple of complications
** introduced by the redundant storage I use for open-shell orbitals
** (useful for spin-restricted references --- see the CCSORT code). In
** particular, if the Mulliken-ordered density is not bra-ket
** symmetric, specific elements of the final density must be
** multiplied by two or they will not appear with the correct
** prefactor in the backtransformation.  This only affects the IJKA,
** IAJB, and ABCI Dirac-ordered densities, since the remaining three
** components are bra-ket symmetric in Mulliken order.
**
** I really need to give an example of this problem using specific
** elements of GIJKA so that the code below will be clearer.*/

void dump_ROHF(struct iwlbuf *OutBuf, struct RHO_Params rho_params)
{
  int nirreps, nmo, nfzv;
  int *qt_occ, *qt_vir;
  int h, row, col, p, q, r, s;
  dpdbuf4 G;

  qt_occ = moinfo.qt_occ;  qt_vir = moinfo.qt_vir;
  nirreps = moinfo.nirreps;
  nmo = moinfo.nmo;
  nfzv = moinfo.nfzv;

  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
 /*  psio_write_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) moinfo.opdm[0], */
  psio_write_entry(PSIF_MO_OPDM, rho_params.opdm_lbl, (char *) moinfo.opdm[0],
		   sizeof(double)*(nmo-nfzv)*(nmo-nfzv));
  psio_close(PSIF_MO_OPDM, 1);

if (!params.onepdm) {
  psio_open(PSIF_MO_LAG, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) moinfo.I[0],
		   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_LAG, 1);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 0, "G(IK,JL)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 0, 0, 0, 0, "G(IK,JL)");
  dpd_buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_occ, 1, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 10, "G(IK,JA)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 10, 0, 10, 0, "G(IK,JA)");
  
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0];
      q = G.params->roworb[h][row][1];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0];
	s = G.params->colorb[h][col][1];

	if((qt_occ[q] == qt_vir[s]) && (p == r))
	  G.matrix[h][row][col] *= 2;
      }
    }

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_vir, 0, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP9, prqs, 10, 10, "G(IA,JB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP9, 0, 10, 10, 10, 10, 0, "G(IA,JB)");
  dpd_buf4_symm(&G);
  dpd_buf4_dump(&G, OutBuf, qt_occ, qt_vir, qt_occ, qt_vir, 1, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 5, "G(IJ,AB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  dpd_buf4_scm(&G, 0.5);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0];
      q = G.params->roworb[h][row][1];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0];
	s = G.params->colorb[h][col][1];

	if((qt_occ[p] == qt_vir[r]) && (qt_occ[q] == qt_vir[s]))
	  G.matrix[h][row][col] *= 2;
      }
    }

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  
  dpd_buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_vir, qt_vir, 0, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 5, 10, "G(ca,IB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 5, 10, 5, 10, 0, "G(ca,IB)");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0];
      q = G.params->roworb[h][row][1];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0];
	s = G.params->colorb[h][col][1];

	if((qt_vir[p] == qt_occ[r]) && (q == s))
	  G.matrix[h][row][col] *= 2;
      }
    }

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_dump(&G, OutBuf, qt_vir, qt_vir, qt_occ, qt_vir, 0, 0);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 5, 5, "G(AC,BD)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 5, 5, 5, 5, 0, "G(AC,BD)");
  dpd_buf4_dump(&G, OutBuf, qt_vir, qt_vir, qt_vir, qt_vir, 1, 0);
  dpd_buf4_close(&G);

  }
}

}} // namespace psi::ccdensity
