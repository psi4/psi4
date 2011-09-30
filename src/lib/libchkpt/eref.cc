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

double Chkpt::rd_eref(void)
{
        double eref;
        char *keyword;
        keyword = build_keyword("Reference energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &eref, sizeof(double));

        free(keyword);
        return eref;
}

void Chkpt::wt_eref(double eref)
{
        char *keyword;
        keyword = build_keyword("Reference energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &eref,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_eref(): Reads in the reference energy.
**
**   takes no arguments.
**
**   returns: double eref  the reference energy.
**
** \ingroup CHKPT
*/
        double chkpt_rd_eref(void)
        {
                double eref;
                eref = _default_chkpt_lib_->rd_eref();
                return eref;
        }


/*!
** chkpt_wt_eref(): Writes out the reference energy.
**
** \param double eref = the reference energy.
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_eref(double eref)
        {
                _default_chkpt_lib_->wt_eref(eref);
        }
}
