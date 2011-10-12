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

double *Chkpt::rd_fock(void)
{
        double *fmat;
        char *keyword;
        keyword = build_keyword("Fock Matrix");

        int nso = rd_nso();

        fmat = array<double>((nso*nso+nso)/2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) fmat,
                           ((nso*nso+nso)/2)*sizeof(double));
        free(keyword);
        return fmat;
}

void Chkpt::wt_fock(double *fmat)
{
        char *keyword;
        keyword = build_keyword("Fock Matrix");

        int nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) fmat,
                ((nso*nso+nso)/2)*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_fock():  Reads in the Fock Matrix
**
**  takes no arguments.
**
**  returns: double *fmat  an array lower triangle closed shell fock matrix
**      ordered by irrep.
** \ingroup CHKPT
*/
        double *chkpt_rd_fock(void)
        {
                double *fmat;
                fmat = _default_chkpt_lib_->rd_fock();
                return fmat;
        }

/*!
** chkpt_wt_fock():  Writes the Fock Matrix
**
** arguments:
**  \param evals = an array of the lower triangle part of the fock matrix
**      ordered by irrep.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fock(double *fmat)
        {
                _default_chkpt_lib_->wt_fock(fmat);
        }

}
