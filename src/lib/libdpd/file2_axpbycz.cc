/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

/* file2_axpbycz(): Evaluates the standard operation aX + bY -> cZ for
** dpdfile2's.
**
** Arguments:
**   dpdfile2 *FileA: A pointer to the leftmost dpdfile2.
**   dpdfile2 *FileB: A pointer to the rightmost summand dpdfile2.
**   dpdfile2 *FileC: A pointer to the target dpdfile2.
**   double a, b, c, scalar prefactors
*/

int dpd_file2_axpbycz(dpdfile2 *FileA, dpdfile2 *FileB, dpdfile2 *FileC,
  double a, double b, double c)
{
  dpd_file2_scm(FileC, c);

  dpd_file2_axpy(FileB, FileC, b, 0);

  dpd_file2_axpy(FileA, FileC, a, 0);
  return 0;
}


} /* extern "C" */
