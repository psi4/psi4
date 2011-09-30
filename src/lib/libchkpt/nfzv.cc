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

int Chkpt::rd_nfzv(void)
{
        int nfzv;
        char *keyword;
        keyword = build_keyword("Num. Frozen UOCC");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int) );

        free(keyword);
        return nfzv;
}

void Chkpt::wt_nfzv(int nfzv)
{
        char *keyword;
        keyword = build_keyword("Num. Frozen UOCC");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nfzv()
** Reads in the total number of frozen unoccupied molecular orbitals.
**
** returns: nfzv = total number of frozen unoccupied molecular orbitals.
** \ingroup CHKPT
*/
        int chkpt_rd_nfzv(void)
        {
                return _default_chkpt_lib_->rd_nfzv();
        }

/*!
** void chkpt_wt_nfzv(int)
** Writes out the total number of frozen unoccupied molecular orbitals.
**
** \param nfzv = total number of frozen unoccupied molecular orbitals.
**
** \ingroup CHKPT
*/
        void chkpt_wt_nfzv(int nfzv)
        {
                _default_chkpt_lib_->wt_nfzv(nfzv);
        }
}
