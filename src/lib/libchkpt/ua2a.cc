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

int *Chkpt::rd_ua2a(void)
{
        int *ua2a;
        int num_unique_atoms;
        char *keyword;
        keyword = build_keyword("Unique atom -> full atom map");

        num_unique_atoms = rd_num_unique_atom();
        ua2a = array<int>(num_unique_atoms);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

        free(keyword);
        return ua2a;
}

void Chkpt::wt_ua2a(int *ua2a)
{
        int num_unique_atoms;
        char *keyword;
        keyword = build_keyword("Unique atom -> full atom map");

        num_unique_atoms = rd_num_unique_atom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_ua2a()
** Reads in a mapping array from the symmetry-unique atom
** list to the full atom list
**
** returns: ua2a = Read in an array num_unique_atom long
**
** \ingroup CHKPT
*/
        int *chkpt_rd_ua2a(void)
        {
                return _default_chkpt_lib_->rd_ua2a();
        }

/*!
** void chkpt_wt_ua2a(int *)
** Writes out a mapping array from the symmetry-unique atom
** list to the full atom list
**
** \param ua2a = An array num_unique_atom long
**
** returns: none
*/
        void chkpt_wt_ua2a(int *ua2a)
        {
                _default_chkpt_lib_->wt_ua2a(ua2a);
        }
}
