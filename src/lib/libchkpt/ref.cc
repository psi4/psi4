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

int Chkpt::rd_ref(void)
{
        int refnum;
        char *keyword;
        keyword = build_keyword("Reference");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &refnum, sizeof(int) );

        free(keyword);
        return refnum;
}

void Chkpt::wt_ref(int refnum)
{
        char *keyword;
        keyword = build_keyword("Reference");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &refnum, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_ref()
** Reads the reference type from the flag in checkpoint
** 0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF
**
** returns: refnum = number indicating the reference.
** \ingroup CHKPT
*/
        int chkpt_rd_ref(void)
        {
                return _default_chkpt_lib_->rd_ref();
        }

/*!
** void chkpt_wt_ref(int)
** Writes out the reference type from the flag in checkpoint
** 0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF
**
** \param refnum = number indicating the reference.
** \ingroup CHKPT
*/
        void chkpt_wt_ref(int refnum)
        {
                _default_chkpt_lib_->wt_ref(refnum);
        }
}
