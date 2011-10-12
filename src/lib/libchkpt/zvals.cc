/*!
    \file
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double *Chkpt::rd_zvals(void)
{
        int natom;
        double *zvals;
        char *keyword;
        keyword = build_keyword("Nuclear charges");

        natom = rd_natom();
        zvals = array<double>(natom);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) zvals,
                natom*sizeof(double));

        free(keyword);
        return zvals;
}

void Chkpt::wt_zvals(double *zvals)
{
        int natom;
        char *keyword;
        keyword = build_keyword("Nuclear charges");

        natom = rd_natom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) zvals,
                natom*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_zvals()
** Reads the nuclear charges from the checkpoint file.
**
** arguments: none
**
** returns:
**   double *zvals: An array of the charges
*/
        double *chkpt_rd_zvals(void)
        {
                return _default_chkpt_lib_->rd_zvals();
        }

/*!
** chkpt_wt_zvals()
** Writes the nuclear charges to the checkpoint file.
**
** \param double *zvals: An array of the charges
**
** returns: nothing
*/
        void chkpt_wt_zvals(double *zvals)
        {
                _default_chkpt_lib_->wt_zvals(zvals);
        }
}
