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

double Chkpt::rd_ecorr(void)
{
        double ecorr;
        char *keyword;
        keyword = build_keyword("Correlation energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
                sizeof(double));

        free(keyword);
        return ecorr;
}

void Chkpt::wt_ecorr(double ecorr)
{
        char *keyword;
        keyword = build_keyword("Correlation energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_ecorr():  Reads in the correlated energy.
**
** takes no arguments.
**
** returns: e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
** \ingroup CHKPT
*/
        double chkpt_rd_ecorr(void)
        {
                double energy;
                energy = _default_chkpt_lib_->rd_ecorr();
                return energy;
        }

/*!
** chkpt_wt_ecorr():  Writes out the correlated energy.
**
** \param e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_ecorr(double ecorr)
        {
                _default_chkpt_lib_->wt_ecorr(ecorr);
        }
}
