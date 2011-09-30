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

int Chkpt::rd_nshell(const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("Num. shells", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nshell, sizeof(int));

        free(keyword);
        return nshell;
}

void Chkpt::wt_nshell(int nshell, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Num. shells", key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nshell, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nshell()
** Reads in the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells.
** Shells of all atoms are counted (compare nprim).
**
** returns: nshell = total number of shells.
** \ingroup CHKPT
*/
        int chkpt_rd_nshell(void)
        {
                return _default_chkpt_lib_->rd_nshell();
        }


/*!
** void chkpt_wt_nshell(int)
** Writes out the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells.
** Shells of all atoms are counted (compare nprim).
**
** \param nshell = total number of shells.
**
** returns:none
**
** \ingroup CHKPT
*/
        void chkpt_wt_nshell(int nshell, const char *key2)
        {
                _default_chkpt_lib_->wt_nshell(nshell, key2);
        }
}
