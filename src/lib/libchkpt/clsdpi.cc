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

int *Chkpt::rd_clsdpi(void)
{
        int nirreps;
        int *clsdpi;
        char *keyword;
        keyword = build_keyword("Closed shells per irrep");

        nirreps = rd_nirreps();
        clsdpi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) clsdpi,
                nirreps*sizeof(int));

        free(keyword);
        return clsdpi;
}

void Chkpt::wt_clsdpi(int *clsdpi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Closed shells per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) clsdpi,
                nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_clsdpi():  Reads in the number of closed-shell orbitals in each
**   irrep.
**
**  Parameters: none
**
**  Returns:
**    int *clsdpi, an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of closed-shell orbitals for
**                 that irrep.
** \ingroup CHKPT
*/
        int *chkpt_rd_clsdpi(void)
        {
                int *clsdpi;
                clsdpi = _default_chkpt_lib_->rd_clsdpi();
                return clsdpi;
        }


/*!
** chkpt_wt_clsdpi():  Writes the number of closed-shell orbitals in each irrep.
**
** \param clsdpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of closed-shell orbitals for
**                 that irrep.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_clsdpi(int *clsdpi)
        {
                _default_chkpt_lib_->wt_clsdpi(clsdpi);
        }
}

