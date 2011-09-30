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

int *Chkpt::rd_openpi(void)
{
        int nirreps;
        int *openpi;
        char *keyword;
        keyword = build_keyword("Open shells per irrep");

        nirreps = rd_nirreps();
        openpi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) openpi,
                nirreps*sizeof(int));

        free(keyword);
        return openpi;
}

void Chkpt::wt_openpi(int *openpi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Open shells per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) openpi,
                nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_openpi(): Reads in the number of open-shell orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     *openpi  an array which has an element for each irrep of the
**              point group of the molecule (n.b. not just the ones
**              with a non-zero number of basis functions). each
**              element contains the number of open-shell orbitals for
**              that irrep.
*/
        int *chkpt_rd_openpi(void)
        {
                return _default_chkpt_lib_->rd_openpi();
        }

/*!
** chkpt_wt_openpi():  Writes the number of open-shell orbitals in each irrep.
**
** \param *openpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of open-shell orbitals for
**                 that irrep.
**
** returns: none
*/
        void chkpt_wt_openpi(int *openpi)
        {
                _default_chkpt_lib_->wt_openpi(openpi);
        }
}
