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

int *Chkpt::rd_frzcpi(void)
{
        int nirreps;
        int *frzcpi;
        char *keyword;
        keyword = build_keyword("Frozen DOCC per irrep");

        nirreps = rd_nirreps();
        frzcpi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) frzcpi,
          nirreps*sizeof(int));

        free(keyword);
        return frzcpi;
}

void Chkpt::wt_frzcpi(int *frzcpi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Frozen DOCC per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) frzcpi,
          nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_frzcpi():  Reads in the number of frozen doubly occupied molecular
**   orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *frzcpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
** \ingroup CHKPT
*/
        int *chkpt_rd_frzcpi(void)
        {
                int *frzcpi;
                frzcpi = _default_chkpt_lib_->rd_frzcpi();
                return frzcpi;
        }


/*!
** chkpt_wt_frzcpi():  Writes the number of frozen doubly occupied molecular
**   orbitals in each irrep
**
** \param frzcpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for that irrep.  See also
**                 chkpt_rd_sopi().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_frzcpi(int *frzcpi)
        {
                _default_chkpt_lib_->wt_frzcpi(frzcpi);
        }
}
