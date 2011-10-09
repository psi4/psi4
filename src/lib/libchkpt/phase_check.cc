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

int Chkpt::rd_phase_check(void)
{
        int pcheck;
        char *keyword;
        keyword = build_keyword("Phase check");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &pcheck, sizeof(int));

        free(keyword);
        return pcheck;
}

void Chkpt::wt_phase_check(int pcheck)
{
        char *keyword;
        keyword = build_keyword("Phase check");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &pcheck, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_phase_check()
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** arguments: none
**
** returns: pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** \ingroup CHKPT
*/
        int chkpt_rd_phase_check(void)
        {
                return _default_chkpt_lib_->rd_phase_check();
        }

/*!
** void chkpt_wt_phase_check(int)
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** \param pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_phase_check(int pcheck)
        {
                _default_chkpt_lib_->wt_phase_check(pcheck);
        }
}
