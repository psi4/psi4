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

int Chkpt::rd_nfzc(void)
{
        int nfzc;
        char *keyword;
        keyword = build_keyword("Num. Frozen DOCC");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nfzc, sizeof(int) );

        free(keyword);
        return nfzc;
}

void Chkpt::wt_nfzc(int nfzc)
{
char *keyword;
keyword = build_keyword("Num. Frozen DOCC");

psio->write_entry(PSIF_CHKPT, keyword, (char *) &nfzc, sizeof(int));

free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nfzc()
** Reads in the total number of frozen doubly occupied molecular orbitals.
**
** returns: nfzc = total number of frozen doubly occupied molecular orbitals.
** \ingroup CHKPT
*/
        int chkpt_rd_nfzc(void)
        {
                return _default_chkpt_lib_->rd_nfzc();
        }

/*!
** void chkpt_wt_nfzc(int)
** Writes out the total number of frozen doubly occupied molecular orbitals.
**
** \param nfzc = total number of frozen doubly occupied molecular orbitals.
**
** \ingroup CHKPT
*/
        void chkpt_wt_nfzc(int nfzc)
        {
                _default_chkpt_lib_->wt_nfzc(nfzc);
        }
}
