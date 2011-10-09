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

int Chkpt::rd_nmo(void)
{
        int nmo;
        char *keyword;
        keyword = build_keyword("Num. MO's");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nmo, sizeof(int));

        free(keyword);
        return nmo;
}

void Chkpt::wt_nmo(int nmo)
{
        char *keyword;
        keyword = build_keyword("Num. MO's");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nmo, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nmo()
** Reads in the total number of molecular orbitals.
**
** returns: nmo = total number of molecular orbitals.
** \ingroup CHKPT
*/
        int chkpt_rd_nmo(void)
        {
                return _default_chkpt_lib_->rd_nmo();
        }

/*!
** void chkpt_wt_nmo(int)
** Writes out the total number of molecular orbitals.
**
** \param nmo = total number of molecular orbitals.
**
** \ingroup CHKPT
*/
        void chkpt_wt_nmo(int nmo)
        {
                _default_chkpt_lib_->wt_nmo(nmo);
        }
}
