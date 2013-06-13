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

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cctriples {

// If gradients are needed, the abc driven loops will be called at some point
// They expect some transposed matrices for efficiency, so here they are.
// ACS June 08
void transpose_integrals() {
  dpdbuf4 T2AA;
  dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_->buf4_sort(&T2AA, PSIF_CC_TAMPS, rspq, 7, 2, "tABIJ");
  dpd_->buf4_close(&T2AA);

  dpdbuf4 T2BB;
  dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  dpd_->buf4_sort(&T2BB, PSIF_CC_TAMPS, rspq, 17, 12, "tabij");
  dpd_->buf4_close(&T2BB);

  dpdbuf4 T2AB;
  dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_->buf4_sort(&T2AB, PSIF_CC_TAMPS, rspq, 28, 22, "tAbIj");
  dpd_->buf4_close(&T2AB);

  dpdbuf4 T2BA;
  dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_->buf4_sort(&T2BA, PSIF_CC_TAMPS, rspq, 29, 23, "taBiJ");
  dpd_->buf4_close(&T2BA);

  dpdbuf4 FAAints;
  dpd_->buf4_init(&FAAints, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_->buf4_sort(&FAAints, PSIF_CC_FINTS, rspq, 7, 20, "F <BC||IA>");
  dpd_->buf4_close(&FAAints);

  dpdbuf4 FBBints;
  dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  dpd_->buf4_sort(&FBBints, PSIF_CC_FINTS, rspq, 17, 30, "F <bc||ia>");
  dpd_->buf4_close(&FBBints);

  dpdbuf4 FABints;
  dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_->buf4_sort(&FABints, PSIF_CC_FINTS, rspq, 28, 24, "F <Bc|Ia>");
  dpd_->buf4_close(&FABints);

  dpdbuf4 FBAints;
  dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  dpd_->buf4_sort(&FBAints, PSIF_CC_FINTS, rspq, 29, 27, "F <bC|iA>");
  dpd_->buf4_close(&FBAints);

  dpdbuf4 EAAints;
  dpd_->buf4_init(&EAAints, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
  dpd_->buf4_sort(&EAAints, PSIF_CC_EINTS, srpq, 21, 2, "E <AK||IJ> (AK, I>J)");
  dpd_->buf4_close(&EAAints);

  dpdbuf4 EBBints;
  dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 11, 30, 11, 30, 0, "E <ij||ka> (i>j,ka)");
  dpd_->buf4_sort(&EBBints, PSIF_CC_EINTS, srpq, 31, 12, "E <ak||ij> (ak, i>j)");
  dpd_->buf4_close(&EBBints);

  dpdbuf4 EABints;
  dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_->buf4_sort(&EABints, PSIF_CC_EINTS, srpq, 25, 22, "E <aK|Ij>");
  dpd_->buf4_close(&EABints);

  dpdbuf4 EBAints;
  dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
  dpd_->buf4_sort(&EBAints, PSIF_CC_EINTS, srpq, 26, 23, "E <Ak|iJ>");
  dpd_->buf4_close(&EBAints);
}

}
}