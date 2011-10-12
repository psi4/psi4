/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_iopen(void)
{
        int iopen;
        char *keyword;
        keyword = build_keyword("Iopen");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &iopen, sizeof(int));

        free(keyword);
        return iopen;
}

void Chkpt::wt_iopen(int iopen)
{
        char *keyword;
        keyword = build_keyword("Iopen");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &iopen, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_iopen()
** Reads in dimensionality of ALPHA and BETA vectors of two-electron
** coupling coefficients for open shells.
**
** Note : IOPEN = MM * (MM + 1), where MM is the total number
** of irreps containing singly occupied orbitals.
**
** returns:
**   iopen = dimensionality of ALPHA and BETA vectors of coupling
**           coefficients for open shells.
** \ingroup CHKPT
*/
        int chkpt_rd_iopen(void)
        {
                return _default_chkpt_lib_->rd_iopen();
        }


/*!
** void chkpt_wt_iopen(int)
** Writes out the dimensionality of ALPHA and BETA vectors of two-electron
** coupling coefficients for open shells.
**
** Note : IOPEN = MM * (MM + 1), where MM is the total number
** of irreps containing singly occupied orbitals.
**
**  arguments:
**   \param iopen = dimensionality of ALPHA and BETA vectors of coupling
**                  coefficients for open shells.
** \ingroup CHKPT
*/
        void chkpt_wt_iopen(int iopen)
        {
                _default_chkpt_lib_->wt_iopen(iopen);
        }
}

