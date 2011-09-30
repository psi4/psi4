/*!
  \file
  \ingroup CHKPT
*/

#include <stdlib.h>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_nfragment(void)
{
        int nfragment;
        char *keyword;
        keyword = build_keyword("Num. fragments");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nfragment, sizeof(int));

        free(keyword);
        return nfragment;
}

void Chkpt::wt_nfragment(int nfragment)
{
        char *keyword;
        keyword = build_keyword("Num. fragments");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nfragment, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nfragment()
** Reads in the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** returns: nfragment = total number of irreducible representations.
** \ingroup CHKPT
*/
        int chkpt_rd_nfragment(void)
        {
                return _default_chkpt_lib_->rd_nfragment();
        }

/*!
** void chkpt_wt_nfragment(int)
** Writes out the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** \param nfragment = total number of irreducible representations.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_nfragment(int nfragment)
        {
                _default_chkpt_lib_->wt_nfragment(nfragment);
        }
}
