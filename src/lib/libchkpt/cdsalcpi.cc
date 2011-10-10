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

int *Chkpt::rd_cdsalcpi(void)
{
        const int nirreps = rd_nirreps();
        int *cdsalcpi = array<int>(nirreps);
        psio_address ptr = PSIO_ZERO;
        char *keyword = build_keyword("cartdisp SALCs per irrep");

        psio->read(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
        return cdsalcpi;
}

void Chkpt::wt_cdsalcpi(const int *cdsalcpi)
{
        const int nirreps = rd_nirreps();
        psio_address ptr = PSIO_ZERO;
        char *keyword = build_keyword("cartdisp SALCs per irrep");

        psio->write(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_cdsalcpi(): Read in number of SALCs per irrep
**
** Parameters: none
**
** Returns: cdsalcpi = An array of nirreps integers.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_cdsalcpi(void)
        {
                return _default_chkpt_lib_->rd_cdsalcpi();
        }

/*!
** chkpt_wt_cdsalcpi(): Writes out number of SALCs per irrep
**
** \param cdsalcpi = An array of nirreps integers
**
** Returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_cdsalcpi(const int *cdsalcpi)
        {
                _default_chkpt_lib_->wt_cdsalcpi(cdsalcpi);
        }
}

