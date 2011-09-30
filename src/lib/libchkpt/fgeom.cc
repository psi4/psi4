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

double **Chkpt::rd_fgeom(void)
{
        int nallatom;
        double **fgeom;
        char *keyword;
        keyword = build_keyword("Full cartesian geometry");

        nallatom = rd_nallatom();
        fgeom = matrix<double>(nallatom,3);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) fgeom[0],
                                        (int) 3*nallatom*sizeof(double));

        free(keyword);
        return  fgeom;
}

void Chkpt::wt_fgeom(double **fgeom)
{
        int nallatom;
        char *keyword;
        keyword = build_keyword("Full cartesian geometry");

        nallatom = rd_nallatom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) fgeom[0],
                (int) 3*nallatom*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_fgeom():  Reads in full cartesian geometry including dummy atoms
**
** takes no arguments.
** returns: double **full_geom;
** \ingroup CHKPT
*/
        double **chkpt_rd_fgeom(void)
        {
                return _default_chkpt_lib_->rd_fgeom();
        }

/*!
** chkpt_wt_fgeom():  Writes out full cartesian geometry including dummy atoms
**
** arguments:
**   \param full_geom = Matrix for cartesian coordinates
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fgeom(double **fgeom)
        {
                _default_chkpt_lib_->wt_fgeom(fgeom);
        }
}
