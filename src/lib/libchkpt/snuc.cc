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

int *Chkpt::rd_snuc(const char *key2)
{
        int *snuc;
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell nucleus", key2);

        nshell = rd_nshell(key2);

        snuc = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

        free(keyword);
        return snuc;
}

void Chkpt::wt_snuc(int *snuc, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell nucleus", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_snuc(): Reads in array of the nuclei numbers shells belong to.
**
**  takes no arguments.
**
**  returns: snuc = an array of the nuclei numbers to which shells
**                  belong to.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_snuc(void)
        {
                return _default_chkpt_lib_->rd_snuc();
        }

/*!
** chkpt_wt_snuc(): Writes out array of the nuclei numbers shells belong to.
**
**  \param snuc = an array of the nuclei numbers to which shells belong to
**
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_snuc(int *snuc, const char *key2)
        {
                _default_chkpt_lib_->wt_snuc(snuc, key2);
        }
}
