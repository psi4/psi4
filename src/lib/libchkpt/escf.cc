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

double Chkpt::rd_escf(void)
{
        double escf;
        char *keyword;
        keyword = build_keyword("SCF energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

        free(keyword);
        return escf;
}

void Chkpt::wt_escf(double escf)
{
        char *keyword;
        keyword = build_keyword("SCF energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &escf, sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_escf(): Reads in the scf energy.
**
**  takes no arguments.
**
**  returns: double escf  the scf energy.
** \ingroup CHKPT
*/
        double chkpt_rd_escf(void)
        {
                double escf;
                escf = _default_chkpt_lib_->rd_escf();
                return escf;
        }

/*!
** chkpt_wt_escf(): Writes out the scf energy.
**
**  arguments:
**   \param double escf  the scf energy.
**
** returns: none
** \ingroup CHKPT
*/

        void chkpt_wt_escf(double escf)
        {
                _default_chkpt_lib_->wt_escf(escf);
        }
}
