/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double Chkpt::rd_efzc(void)
{
        double efzc;
        char *keyword;
        keyword = build_keyword("Frozen core energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &efzc, sizeof(double));

        free(keyword);
        return efzc;
}

void Chkpt::wt_efzc(double efzc)
{
        char *keyword;
        keyword = build_keyword("Frozen core energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &efzc,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_efzc(): Reads in the frozen-core energy.
**
**   takes no arguments.
**
**   returns: double efzc  the frozen-core energy.
** \ingroup CHKPT
*/
        double chkpt_rd_efzc(void)
        {
                double efzc;
                efzc = _default_chkpt_lib_->rd_efzc();
                return efzc;
        }

/*!
** chkpt_wt_efzc(): Writes out the frozen-core energy.
**
** \param efzc = the frozen-core energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_efzc(double efzc)
        {
                _default_chkpt_lib_->wt_efzc(efzc);
        }
}
