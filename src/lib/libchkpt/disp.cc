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

int Chkpt::rd_disp(void)
{
        int disp;
        char *keyword;
        keyword = build_keyword("Current displacement");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &disp,
                sizeof(int));

        free(keyword);
        return disp;
}

void Chkpt::wt_disp(int disp)
{
        char *keyword;
        keyword = build_keyword("Current displacement");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &disp,
                sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_disp():  Reads in the current geometry displacement number.
**
** Parameters: none
**
** not used by OPTKING; used by anybody else ???
**
** Returns: int disp, the current geometry displacement number
** \ingroup CHKPT
*/
        int chkpt_rd_disp(void)
        {
                int disp;
                disp = _default_chkpt_lib_->rd_disp();
                return disp;
        }


/*!
** chkpt_wt_disp():  Writes out the current geometry displacement number.
**
** \param int disp   the current geometry displacement number
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_disp(int disp)
        {
                _default_chkpt_lib_->wt_disp(disp);
        }
}

