/*!
    \file Read and write 3 rotational constants
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double *Chkpt::rd_rotconst(void)
{
  double *rotconst;
  char *keyword;
  keyword = build_keyword("Rotational Constants");

  rotconst = array<double>(3);

  psio->read_entry(PSIF_CHKPT, keyword, (char *) rotconst,
    3*sizeof(double));

  free(keyword);
  return rotconst;
}

void Chkpt::wt_rotconst(double *rotconst)
{
  char *keyword;
  keyword = build_keyword("Rotational Constants");

  psio->write_entry(PSIF_CHKPT, keyword, (char *) rotconst,
    3*sizeof(double));

  free(keyword);
}

extern "C" {
/*!
** chkpt_rd_rotconst()
** Reads the rotational constants from the checkpoint file.
**
** arguments: none
**
** returns:
**   double *rotconst: An array of the rotational constants
*/
  double *chkpt_rd_rotconst(void)
  {
    return _default_chkpt_lib_->rd_rotconst();
  }

/*!
** chkpt_wt_rotconst()
** Writes the rotational constants to the checkpoint file.
**
** \param double *rotconst: An array of the rotational constants
**
** returns: nothing
*/
  void chkpt_wt_rotconst(double *rotconst)
  {
    _default_chkpt_lib_->wt_rotconst(rotconst);
  }
}
