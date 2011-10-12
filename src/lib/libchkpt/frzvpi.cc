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

int *Chkpt::rd_frzvpi(void)
{
        int nirreps;
        int *frzvpi;
        char *keyword;
        keyword = build_keyword("Frozen UOCC per irrep");

        nirreps = rd_nirreps();
        frzvpi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) frzvpi,
                nirreps*sizeof(int));

        free(keyword);
        return frzvpi;
}

void Chkpt::wt_frzvpi(int *frzvpi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Frozen UOCC per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) frzvpi,
          nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_frzvpi():  Reads in the number of frozen unoccupied molecular
**   orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *frzvpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen unoccupied
**                 molecular orbitals for that irrep.
**                 See also chkpt_rd_sopi().
** \ingroup CHKPT
*/
        int *chkpt_rd_frzvpi(void)
        {
                return _default_chkpt_lib_->rd_frzvpi();
        }

/*!
** chkpt_wt_frzvpi():  Writes the number of frozen unoccupied molecular
**   orbitals in each irrep.
**
** \param frzvpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions).  Each element
**                 contains the number of frozen unoccupied molecular orbitals
**                 for that irrep.   See also chkpt_rd_sopi().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_frzvpi(int *frzvpi)
        {
                _default_chkpt_lib_->wt_frzvpi(frzvpi);
        }
}

