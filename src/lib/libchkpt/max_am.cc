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

int Chkpt::rd_max_am(const char *key2)
{
        int max_am;
        char *keyword;
        keyword = build_keyword("Max. AM", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

        free(keyword);
        return max_am;
}

void Chkpt::wt_max_am(int max_am, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Max. AM", key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_max_am()
** Reads in the maximum orbital quantum number
** of AOs in the basis.
**
** Returns: max_am = (0 corresponds to s-functions,
**                    1 - to up to p-functions, etc.)
** \ingroup CHKPT
*/
        int chkpt_rd_max_am(void)
        {
                int max_am;
                max_am = _default_chkpt_lib_->rd_max_am();
                return max_am;
        }

/*!
** void chkpt_wt_max_am()
** Writes out the maximum orbital quantum number
** of AOs in the basis.
**
** \param max_am = (0 corresponds to s-functions,
**                   1 - to up to p-functions, etc.)
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_max_am(int max_am, const char *key2)
        {
                _default_chkpt_lib_->wt_max_am(max_am, key2);
        }
}

