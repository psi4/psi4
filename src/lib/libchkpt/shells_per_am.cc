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

int *Chkpt::rd_shells_per_am(const char *key2)
{
        int *shells_per_am;
        int max_am;
        char *keyword;
        keyword = build_keyword("Shells per am", key2);

        max_am = rd_max_am(key2);
        shells_per_am = array<int>(max_am+1);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) shells_per_am,
                (max_am+1)*sizeof(int));

        free(keyword);
        return shells_per_am;
}

void Chkpt::wt_shells_per_am(int *shells_per_am, const char *key2)
{
        int max_am;
        char *keyword;
        keyword = build_keyword("Shells per am", key2);

        max_am = rd_max_am(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) shells_per_am,
                (max_am+1)*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_shells_per_am()
** Reads in the numbers of shells of each angular momentum.
**
** returns: shells_per_am = array of shells per angular momentum
**
** \ingroup CHKPT
*/
        int *chkpt_rd_shells_per_am(void)
        {
                return _default_chkpt_lib_->rd_shells_per_am();
        }

/*!
** void chkpt_wt_shells_per_am(int *)
** Writes out the numbers of shells of each angular momentum.
**
** \param shells_per_am = array of shells per angular momentum
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_shells_per_am(int *shells_per_am, const char *key2)
        {
                _default_chkpt_lib_->wt_shells_per_am(shells_per_am,key2);
        }
}
