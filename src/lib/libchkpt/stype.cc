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

int *Chkpt::rd_stype(const char *key2)
{
        int *stype;
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell ang. mom.", key2);

        nshell = rd_nshell(key2);

        stype = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) stype, nshell*sizeof(int));

        free(keyword);
        return stype;
}

void Chkpt::wt_stype(int *stype, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell ang. mom.", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) stype, nshell*sizeof(int));

        free(keyword);
}


extern "C" {
/*!
** chkpt_rd_stype(): 	Reads in an array of the angular momentum numbers of
**			the shells.
**
**  takes no arguments.
**
**  returns: stype = an array of the angular momentum numbers of the shells
**
** \ingroup CHKPT
*/
        int *chkpt_rd_stype(void)
        {
                return _default_chkpt_lib_->rd_stype();
        }

/*!
** chkpt_wt_stype(): 	Writes out an array of the angular momentum numbers of
**			the shells.
**
**  \param stype = an array of the angular momentum numbers of the shells
**
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_stype(int *stype, const char *key2)
        {
                _default_chkpt_lib_->wt_stype(stype, key2);
        }
}
