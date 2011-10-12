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

double *Chkpt::rd_grad(void)
{
        int natom;
        double *grad;
        char *keyword;
        keyword = build_keyword("Energy Gradient");

        natom = rd_natom();
        grad = array<double>(natom*3);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) grad,
          natom*3*sizeof(double));

        free(keyword);
        return grad;
}

void Chkpt::wt_grad(double *grad)
{
        int natom;
        char *keyword;
        keyword = build_keyword("Energy Gradient");

        natom = rd_natom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) grad,
          natom*3*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_grad():  Reads the energy gradient WRT nuclear coordinates
**
**   takes no arguments.
**
**   returns: grad = a vector of doubles natom*3 elements long, e.g.
**     grad[0] = gradient wrt x coordinate of atom 0
**     grad[1] = gradient wrt y coordinate of atom 0
**     grad[8] = gradient wrt z coordinate of atom 2
** \ingroup CHKPT
*/
        double *chkpt_rd_grad(void)
        {
                return _default_chkpt_lib_->rd_grad();
        }

/*!
** chkpt_wt_grad():  Writes the energy gradient WRT nuclear coordinates
**
**   arguments:
**   \param grad = a vector of doubles natom*3 elements long, e.g.
**     grad[0] = gradient wrt x coordinate of atom 0
**     grad[1] = gradient wrt y coordinate of atom 0
**     grad[8] = gradient wrt z coordinate of atom 2
**
**   returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_grad(double *grad)
        {
                _default_chkpt_lib_->wt_grad(grad);
        }
}

