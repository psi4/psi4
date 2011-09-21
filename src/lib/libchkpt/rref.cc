/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double **Chkpt::rd_rref(void)
{
//	char *key;
        double **Rref;
        char *keyword;
        keyword = build_keyword("Transmat to reference frame");

        Rref = matrix<double>(3,3);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

        free(keyword);
        return Rref;
}

void Chkpt::wt_rref(double **Rref)
{
        char *keyword;
        keyword = build_keyword("Transmat to reference frame");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_rref()
** Reads in a 3x3 matrix used to rotate back to the reference frame.
**
**   takes no arguments.
**
**   returns: rref = A 3x3 matrix describing the rotation back to the
**            reference frame.  The reference frame is a coordinate system
**            defined by the "raw" geometry specification (either Z-matrix
**            or geometry array in input.dat or chkpt). Can be used to
**            transform quantities corresponding to different but similar
**            calculations (gradients at displaced geometries) to a
**            common frame.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_rref(void)
        {
                return _default_chkpt_lib_->rd_rref();
        }

/*!
** chkpt_wt_rref()
** Writes out a 3x3 matrix used to rotate back to the reference frame.
**
** \params rref = A 3x3 matrix describing the rotation back to the reference
**                frame.  The reference frame is a coordinate system defined
**                by the "raw" geometry specification (either Z-matrix or
**                geometry array in input.dat or chkpt). Can be used to
**                transform quantities corresponding to different but
**                similar calculations (gradients at displaced geometries)
**                to a common frame.
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_rref(double **Rref)
        {
                _default_chkpt_lib_->wt_rref(Rref);
        }
}
