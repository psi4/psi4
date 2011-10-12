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

int *Chkpt::rd_snumg(const char *key2)
{
        int *snumg;
        int nshell;
        char *keyword;
        keyword = build_keyword("Primitives per shell", key2);

        nshell = rd_nshell(key2);
        snumg = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) snumg, nshell*sizeof(int));

        free(keyword);
        return snumg;
}

void Chkpt::wt_snumg(int *snumg, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("Primitives per shell", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) snumg, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_snumg()
**
** Reads in array of the numbers of the primitive Gaussians in shells.
**
**  takes no arguments.
**
**  returns:
**    snumg = Reads in array of the numbers of the primitive Gaussians
**            in shells
**
** \ingroup CHKPT
*/
        int *chkpt_rd_snumg(void)
        {
                return _default_chkpt_lib_->rd_snumg();
        }

/*!
** chkpt_wt_snumg()
**
** Writes out array of the numbers of the primitive Gaussians in shells.
**
**  \param snumg = array of the numbers of the primitive Gaussians
**                 in shells
**
** \ingroup CHKPT
*/
        void chkpt_wt_snumg(int *snumg, const char *key2)
        {
                _default_chkpt_lib_->wt_snumg(snumg, key2);
        }
}
