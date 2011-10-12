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

int Chkpt::rd_rottype(void)
{
        int rottype;
        char *keyword;
        keyword = build_keyword("Rotor type");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &rottype, sizeof(int));

        free(keyword);
        return rottype;
}

void Chkpt::wt_rottype(int rottype)
{
        char *keyword;
        keyword = build_keyword("Rotor type");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &rottype, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_rottype()
** Reads in type of the rigid rotor molecule represents.
**
** returns: rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
** \ingroup CHKPT
*/
        int chkpt_rd_rottype(void)
        {
                return _default_chkpt_lib_->rd_rottype();
        }

/*!
** void chkpt_wt_rottype(int)
** Reads in type of the rigid rotor molecule represents.
**
** \param rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_rottype(int rottype)
        {
                _default_chkpt_lib_->wt_rottype(rottype);
        }
}
