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

int Chkpt::rd_nsymhf(void)
{
        int nsymhf;
        char *keyword;
        keyword = build_keyword("Num. HF irreps");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

        free(keyword);
        return nsymhf;
}

void Chkpt::wt_nsymhf(int nsymhf)
{
        char *keyword;
        keyword = build_keyword("Num. HF irreps");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nsymhf()
** Reads in the total number of irreps in the point group
** in which the molecule is being considered which
** have non-zero number of basis functions.
**
** returns: nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup CHKPT
*/
        int chkpt_rd_nsymhf(void)
        {
                return _default_chkpt_lib_->rd_nsymhf();
        }

/*!
** void chkpt_wt_nsymhf(int)
** Writes out the total number of irreps in the point group
** in which the molecule is being considered which
** have non-zero number of basis functions.
**
** \param nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup CHKPT
*/
        void chkpt_wt_nsymhf(int nsymhf)
        {
                _default_chkpt_lib_->wt_nsymhf(nsymhf);
        }
}
