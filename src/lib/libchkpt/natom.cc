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

int Chkpt::rd_natom(void)
{
        int natom;
        char *keyword;
        keyword = build_keyword("Num. atoms");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &natom, sizeof(int));

        free(keyword);
        return natom;
}

void Chkpt::wt_natom(int natom)
{
        char *keyword;
        keyword = build_keyword("Num. atoms");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &natom, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_natom()
** Reads in the total number of atoms.
**
** Parameters: none
**
** Returns:
**   natom = total number of atoms.
** \ingroup CHKPT
*/
        int chkpt_rd_natom(void)
        {
                return _default_chkpt_lib_->rd_natom();
        }

/*!
** void chkpt_wt_natom(int)
** Writes out the total number of atoms.
**
** Parameters:
**   \param natom = total number of atoms.
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_natom(int natom)
        {
                _default_chkpt_lib_->wt_natom(natom);
        }
}

