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

int *Chkpt::rd_sprim(const char *key2)
{
        int *sprim;
        int nshell;
        char *keyword;
        keyword = build_keyword("First primitive per shell", key2);

        nshell = rd_nshell(key2);

        sprim = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

        free(keyword);
        return sprim;
}

void Chkpt::wt_sprim(int *sprim, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("First primitive per shell", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_sprim(): Reads in array of the numbers of first primitives
**                   from the shells.
**
**  takes no arguments.
**
**  returns: sprim = an array of the numbers of first primitives
**           from the shells.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_sprim(void)
        {
                return _default_chkpt_lib_->rd_sprim();
        }

/*!
** chkpt_wt_sprim():	Writes out an array of the numbers of first primitives
**			from the shells.
**
**  \param sprim = an array of the numbers of first primitives
**                 from the shells.
**
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_sprim(int *sprim, const char *key2)
        {
                _default_chkpt_lib_->wt_sprim(sprim, key2);
        }
}
