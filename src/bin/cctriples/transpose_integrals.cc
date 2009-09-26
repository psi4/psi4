#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace CCTRIPLES {

// If gradients are needed, the abc driven loops will be called at some point
// They expect some transposed matrices for efficiency, so here they are.
// ACS June 08
void transpose_integrals() {
  dpdbuf4 T2AA;
  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_sort(&T2AA, CC_TAMPS, rspq, 7, 2, "tABIJ");
  dpd_buf4_close(&T2AA);

  dpdbuf4 T2BB;
  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  dpd_buf4_sort(&T2BB, CC_TAMPS, rspq, 17, 12, "tabij");
  dpd_buf4_close(&T2BB);

  dpdbuf4 T2AB;
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_sort(&T2AB, CC_TAMPS, rspq, 28, 22, "tAbIj");
  dpd_buf4_close(&T2AB);

  dpdbuf4 T2BA;
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_buf4_sort(&T2BA, CC_TAMPS, rspq, 29, 23, "taBiJ");
  dpd_buf4_close(&T2BA);

  dpdbuf4 FAAints;
  dpd_buf4_init(&FAAints, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_buf4_sort(&FAAints, CC_FINTS, rspq, 7, 20, "F <BC||IA>");
  dpd_buf4_close(&FAAints);

  dpdbuf4 FBBints;
  dpd_buf4_init(&FBBints, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  dpd_buf4_sort(&FBBints, CC_FINTS, rspq, 17, 30, "F <bc||ia>");
  dpd_buf4_close(&FBBints);

  dpdbuf4 FABints;
  dpd_buf4_init(&FABints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_buf4_sort(&FABints, CC_FINTS, rspq, 28, 24, "F <Bc|Ia>");
  dpd_buf4_close(&FABints);

  dpdbuf4 FBAints;
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  dpd_buf4_sort(&FBAints, CC_FINTS, rspq, 29, 27, "F <bC|iA>");
  dpd_buf4_close(&FBAints);

  dpdbuf4 EAAints;
  dpd_buf4_init(&EAAints, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
  dpd_buf4_sort(&EAAints, CC_EINTS, srpq, 21, 2, "E <AK||IJ> (AK, I>J)");
  dpd_buf4_close(&EAAints);

  dpdbuf4 EBBints;
  dpd_buf4_init(&EBBints, CC_EINTS, 0, 11, 30, 11, 30, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_sort(&EBBints, CC_EINTS, srpq, 31, 12, "E <ak||ij> (ak, i>j)");
  dpd_buf4_close(&EBBints);

  dpdbuf4 EABints;
  dpd_buf4_init(&EABints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_buf4_sort(&EABints, CC_EINTS, srpq, 25, 22, "E <aK|Ij>");
  dpd_buf4_close(&EABints);

  dpdbuf4 EBAints;
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
  dpd_buf4_sort(&EBAints, CC_EINTS, srpq, 26, 23, "E <Ak|iJ>");
  dpd_buf4_close(&EBAints);
}

}
}