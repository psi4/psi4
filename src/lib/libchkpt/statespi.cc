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

int *Chkpt::rd_statespi(void)
{
        int nirreps;
        int *statespi;
        char *keyword;
        keyword = build_keyword("States per irrep");

        nirreps = rd_nirreps();
        statespi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) statespi,
                nirreps*sizeof(int));

        free(keyword);
        return statespi;
}

void Chkpt::wt_statespi(int *statespi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("States per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) statespi,
                nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_statespi():  Reads in the number of excited-states for each
**   irrep.
**
**   takes no arguments.
**
**   returns:
**     int *statespi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of excited states of that
**                 irrep to be studied.
** \ingroup CHKPT
*/
        int *chkpt_rd_statespi(void)
        {
                return _default_chkpt_lib_->rd_statespi();
        }

/*!
** chkpt_wt_statespi():  Writes the number of excited states in each irrep.
**
** \param statespi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of excited states of that
**                 irrep to be studied.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_statespi(int *statespi)
        {
                _default_chkpt_lib_->wt_statespi(statespi);
        }
}
