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

int Chkpt::rd_num_unique_atom(void)
{
        int nunique;
        char *keyword;
        keyword = build_keyword("Num. unique atoms");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

        free(keyword);
        return nunique;
}

void Chkpt::wt_num_unique_atom(int nunique)
{
char *keyword;
keyword = build_keyword("Num. unique atoms");

psio->write_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_num_unique_atom()
** Reads in the number of symmetry unique atoms.
**
** returns: nunique = number of symmetry unique atoms.
** \ingroup CHKPT
*/
        int chkpt_rd_num_unique_atom(void)
        {
                return _default_chkpt_lib_->rd_num_unique_atom();
        }

/*!
** void chkpt_wt_num_unique_atom(int)
** Writes out the number of symmetry unique atoms.
**
** \param nunique = number of symmetry unique atoms.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_num_unique_atom(int nunique)
        {
                _default_chkpt_lib_->wt_num_unique_atom(nunique);
        }
}
