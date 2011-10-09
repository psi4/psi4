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

int *Chkpt::rd_us2s(const char *key2)
{
        int *us2s;
        int num_unique_shells;
        char *keyword;
        keyword = build_keyword("Unique shell -> full shell map", key2);

        num_unique_shells = rd_num_unique_shell(key2);
        us2s = array<int>(num_unique_shells);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

        free(keyword);
        return us2s;
}

void Chkpt::wt_us2s(int *us2s, const char *key2)
{
        int num_unique_shells;
        char *keyword;
        keyword = build_keyword("Unique shell -> full shell map", key2);

        num_unique_shells = rd_num_unique_shell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_us2s()
** Read in a mapping array betwen unique shell and
** full shell lists
**
** returns: us2s = Read in an array num_unique_shell
**
** \ingroup CHKPT
*/
        int *chkpt_rd_us2s(void)
        {
                return _default_chkpt_lib_->rd_us2s();
        }

/*!
** void chkpt_wt_us2s(int *)
** Writes out a mapping array betwen unique shell and
** full shell lists.
**
** \param us2s = An array num_unique_shell
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_us2s(int *us2s, const char *key2)
        {
                _default_chkpt_lib_->wt_us2s(us2s, key2);
        }
}
