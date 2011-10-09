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

double Chkpt::rd_eccsd()
{
        double energy;
        char *keyword;
        keyword = build_keyword("CCSD Energy");

        // Read the energy in
        psio->read_entry(PSIF_CHKPT, keyword, (char*)&energy, sizeof(double));

        // Return the value to the user
        return energy;
}

void Chkpt::wt_eccsd(double eccsd)
{
        char *keyword;
        keyword = build_keyword("CCSD Energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char*)&eccsd, sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_eccsd(): Reads in the CCSD contribution to total energy.
**
**   takes no arguments.
**
**   returns: double eccsd  the CCSD energy.
** \ingroup CHKPT
*/
        double chkpt_rd_eccsd(void)
        {
                double eccsd;
                eccsd = _default_chkpt_lib_->rd_eccsd();
                return eccsd;
        }

/*!
** chkpt_wt_eccsd(): Writes out the CCSD contribution to total energy.
**
** \param eccsd = the CCSD energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_eccsd(double eccsd)
        {
                _default_chkpt_lib_->wt_eccsd(eccsd);
        }
}
