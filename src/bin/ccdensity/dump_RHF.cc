/*! \file 
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
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

    /* DUMP_RHF(): Mulliken-order the RHF-CC two-electron density and
    ** dump it to a file for subsequent backtransformation.  Basically
    ** all we have to do is swap indices two and three, e.g.
    **
    ** G'(pr,qs) = G(pq,rs)
    **
    ** In order for the Mulliken-ordered density to be valid for the
    ** backtransformation algorithm used in TRANSQT, the final density
    ** must have eight-fold permutational symmetry like the original
    ** integrals, so symmetrization of some components is necessary.
    **
    ** Note that the CC-ordered indices must be translated inside
    ** dpd_buf4_dump() to QT ordering, as required by TRANSQT, hence
    ** the qt_occ and qt_vir arrays.
    **
    ** The bk_pack option to dpd_buf4_dump() indicates whether the
    ** Mulliken-ordered twopdm component is bra-ket symmetric
    ** already.  If so, then only the lower-triangle needs to be sent
    ** to TRANSQT2.
    **
    ** I don't remember why I have the swap23 option in
    ** dpd_buf4_dump(), particularly since I do an explicit sort of
    ** the density components using DPD functions.  It would actually
    ** be better to let the dump function handle this part to reduce
    ** disk usage.
    **
    ** Finally, I note that Gibja is multiplied by 0.5 before it is
    ** dumped to disk.  I don't remember why this is necessary.
    **
    ** TDC, last updated 2/08
    */

    void dump_RHF(struct iwlbuf *OutBuf, struct RHO_Params rho_params)
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
	dpd_buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_vir, 0, 0);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	dpd_buf4_sort(&G, CC_TMP9, prqs, 10, 10, "G(IA,JB)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP9, 0, 10, 10, 10, 10, 0, "G(IA,JB)");
	dpd_buf4_symm(&G);
	dpd_buf4_dump(&G, OutBuf, qt_occ, qt_vir, qt_occ, qt_vir, 1, 0);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 5, "G(IJ,AB)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
	//        dpd_buf4_scm(&G, 0.5);  /* why is this necessary? */
	dpd_buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_vir, qt_vir, 0, 0);
	dpd_buf4_close(&G);

	dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	dpd_buf4_sort(&G, CC_TMP0, prqs, 5, 10, "G(ca,IB)");
	dpd_buf4_close(&G);
	dpd_buf4_init(&G, CC_TMP0, 0, 5, 10, 5, 10, 0, "G(ca,IB)");
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
