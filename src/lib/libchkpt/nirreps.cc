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

int Chkpt::rd_nirreps(void)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Num. irreps");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nirreps, sizeof(int));

        free(keyword);
        return nirreps;
}

void Chkpt::wt_nirreps(int nirreps)
{
        char *keyword;
        keyword = build_keyword("Num. irreps");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nirreps, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nirreps()
** Reads in the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** returns: nirreps = total number of irreducible representations.
** \ingroup CHKPT
*/
        int chkpt_rd_nirreps(void)
        {
                return _default_chkpt_lib_->rd_nirreps();
        }

/*!
** void chkpt_wt_nirreps(int)
** Writes out the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** \param nirreps = total number of irreducible representations.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_nirreps(int nirreps)
        {
                _default_chkpt_lib_->wt_nirreps(nirreps);
        }
}
