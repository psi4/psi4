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

int Chkpt::rd_ncalcs(void)
{
        char *keyword_mo, *keyword_alpha_mo;
        keyword_mo = build_keyword("MO coefficients");
        keyword_alpha_mo = build_keyword("Alpha MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword_mo) == NULL &&
                psio->tocscan(PSIF_CHKPT, keyword_alpha_mo) == NULL)
                return 0;
        else
                return 1;

        free(keyword_mo);
        free(keyword_alpha_mo);
}

extern "C" {
/*!
** int chkpt_rd_ncalcs()
** Reads in the total number of HF wave functions.
**
** returns: ncalcs = total number of HF wave functions in checkpoint
** \ingroup CHKPT
*/
        int chkpt_rd_ncalcs(void)
        {
                return _default_chkpt_lib_->rd_ncalcs();
        }
}

