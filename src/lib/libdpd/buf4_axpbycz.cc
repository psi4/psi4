/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {
	
/* buf4_axpbycz(): Evaluates the standard operation aX + bY -> cZ for
** dpdbuf4's.
**
** Arguments:
**   dpdbuf4 *FileA: A pointer to the leftmost dpdbuf4.
**   dpdbuf4 *FileB: A pointer to the rightmost summand dpdbuf4.
**   dpdbuf4 *FileC: A pointer to the target dpdbuf4.
**   double a, b, c, scalar prefactors
*/

int dpd_buf4_axpbycz(dpdbuf4 *FileA, dpdbuf4 *FileB, dpdbuf4 *FileC,
		     double a, double b, double c)
{
  dpd_buf4_scm(FileC, c);

  dpd_buf4_axpy(FileB, FileC, b);

  dpd_buf4_axpy(FileA, FileC, a);
  return 0;
}

} // namespace psi
