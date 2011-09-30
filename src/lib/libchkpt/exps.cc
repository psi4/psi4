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

double *Chkpt::rd_exps(const char *key2)
{
        double *exps;
        int nprim = 0;
        char *keyword;
        keyword = build_keyword("Exponents",key2);

        nprim = rd_nprim(key2);
        exps = array<double>(nprim);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) exps,
                nprim*sizeof(double));

        free(keyword);
        return exps;
}

void Chkpt::wt_exps(double *exps, const char *key2)
{
        int nprim;
        char *keyword;
        keyword = build_keyword("Exponents",key2);

        nprim = rd_nprim(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) exps,
                nprim*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_exps():
** Reads in the exponents of the primitive Gaussian functions.
**
** takes no arguments.
**
** returns: double *exps
** The exponents are returned as an array of doubles.
** \ingroup CHKPT
*/
        double *chkpt_rd_exps(void)
        {
                double *exps;
                exps = _default_chkpt_lib_->rd_exps();
                return exps;
        }

/*!
** chkpt_wt_exps():
** Writes out the exponents of the primitive Gaussian functions.
**
** arguments:
**  \param exps = The exponents are returned as an array of doubles.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_exps(double *exps, const char *key2)
        {
                _default_chkpt_lib_->wt_exps(exps, key2);
        }
}
