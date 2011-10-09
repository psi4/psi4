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

int Chkpt::rd_num_unique_shell(const char *key2)
{
        int nunique;
        char *keyword;
        keyword = build_keyword("Num. unique shells", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

        free(keyword);
        return nunique;
}

void Chkpt::wt_num_unique_shell(int nunique, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Num. unique shells", key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_num_unique_shell()
** Reads in the number of symmetry unique shells.
**
** returns: nunique = number of symmetry unique shells.
** \ingroup CHKPT
*/
        int chkpt_rd_num_unique_shell(void)
        {
                return _default_chkpt_lib_->rd_num_unique_shell();
        }

/*!
** void chkpt_wt_num_unique_shell(int)
** Writes out the number of symmetry unique shells.
**
** \param nunique = number of symmetry unique shells.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_num_unique_shell(int nunique, const char *key2)
        {
                _default_chkpt_lib_->wt_num_unique_shell(nunique, key2);
        }
}
