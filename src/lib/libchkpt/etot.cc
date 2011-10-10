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

double Chkpt::rd_etot(void)
{
        double etot;
        char *keyword;
        keyword = build_keyword("Total energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

        free(keyword);
        return etot;
}

void Chkpt::wt_etot(double etot)
{
        char *keyword;
        keyword = build_keyword("Total energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_etot(): Reads in the total energy.
**
**  takes no arguments.
**
**  returns: double etot  the total energy.
**  \ingroup CHKPT
*/
        double chkpt_rd_etot(void)
        {
                double etot;
                etot = _default_chkpt_lib_->rd_etot();
                return etot;
        }

/*!
** chkpt_wt_etot(): Writes out the total energy.
**
**  arguments:
**   \param double etot  the total energy.
**
**  returns: none
**  \ingroup CHKPT
*/
        void chkpt_wt_etot(double etot)
        {
                _default_chkpt_lib_->wt_etot(etot);
        }
}
