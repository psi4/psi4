/*!
  \file
  \ingroup CHKPT
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <libciomr/libciomr.h>

using namespace psi;

double ***Chkpt::rd_fragment_coeff(void)
{
        int nfragment, *natom_per_fragment, *nref_per_fragment, i, j;
    psio_address ptr;
        double ***fragment_coeff;
        char *keyword;
        keyword = build_keyword("Fragment coeff");

        nfragment = rd_nfragment();
        natom_per_fragment = rd_natom_per_fragment();
        nref_per_fragment = rd_nref_per_fragment();

        fragment_coeff = array<double **>(nfragment);
    for (i=0; i<nfragment; ++i)
          fragment_coeff[i] = matrix<double>(nref_per_fragment[i],natom_per_fragment[i]);

    ptr = PSIO_ZERO;
    for (i=0; i<nfragment; ++i)
      for (j=0; j<nref_per_fragment[i]; ++j)
            psio->read(PSIF_CHKPT, keyword, (char *) fragment_coeff[i][j],
            (int) natom_per_fragment[i]*sizeof(double), ptr, &ptr);

    free(natom_per_fragment);
    free(nref_per_fragment);
        free(keyword);
        return fragment_coeff;
}

void Chkpt::wt_fragment_coeff(double ***fragment_coeff)
{
        int nfragment, *natom_per_fragment, *nref_per_fragment, i, j;
    psio_address ptr;
        char *keyword;
        keyword = build_keyword("Fragment coeff");

        nfragment = rd_nfragment();
    natom_per_fragment = rd_natom_per_fragment();
    nref_per_fragment = rd_nref_per_fragment();

    ptr = PSIO_ZERO;
    for (i=0; i<nfragment; ++i)
      for (j=0; j<nref_per_fragment[i]; ++j)
            psio->write(PSIF_CHKPT, keyword, (char *) fragment_coeff[i][j],
            (int) natom_per_fragment[i]*sizeof(double), ptr, &ptr);

    free(natom_per_fragment);
    free(nref_per_fragment);
        free(keyword);
    return;
}

extern "C" {
/*!
** chkpt_rd_fragment_coeff():  Reads in the coefficients specifying reference points
** for molecular fragments
**
**   takes no arguments.
**
**   returns:
**     double ***fragment_coeff[fragment][reference point][atom in fragment]
** \ingroup CHKPT
*/
        double ***chkpt_rd_fragment_coeff(void)
        {
                double ***fragment_coeff;
                fragment_coeff = _default_chkpt_lib_->rd_fragment_coeff();
                return fragment_coeff;
        }


/*!
** chkpt_wt_fragment_coeff():  Writes out the coefficients specifying the reference
** points for molecular fragments
**
** \param double ***fragment_coeff[fragment][reference point][atom in fragment]
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fragment_coeff(double ***fragment_coeff)
        {
                _default_chkpt_lib_->wt_fragment_coeff(fragment_coeff);
        }
}
