/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

/* build_B_RHF(): Builds the RHF B matrix for RPA calculations.
** In spin orbitals, the B matrix is:
**
** B(ai,bj) = <ij||ab>
**
** RHF references and singlet eigenstates:
**  B(AI,BJ) = 2 <IJ|AB> - <IJ|BA>
**
** RHF references and triplet eigenstates:
**  B(AI,BJ) = <IJ|AB>
**
** TDC, March 2003
*/

void build_B_RHF(void)
{
  dpdbuf4 D;

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "B(AI,BJ)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>"); 
  dpd_buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "B(AI,BJ) triplet");
  dpd_buf4_close(&D);
}

}} // namespace psi::response
