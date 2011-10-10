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

double *Chkpt::rd_contr(const char */*key2*/)
{
#if 0
        double *contr;
        double *temp_contr;
        int nprim, i, j, ij = 0;
        char *keyword;
        keyword = build_keyword("Contraction coefficients", key2);

        nprim = rd_nprim(key2);

        temp_contr = array<double>(MAXANGMOM*nprim);
        contr = array<double>(nprim);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) temp_contr,
                MAXANGMOM*nprim*sizeof(double));

/* Picking non-zero coefficients to the "master" array contr */
        for(i=0; i < MAXANGMOM; i++)
                for(j=0; j < nprim; j++)
        {
                if (temp_contr[ij] != 0)
                        contr[j] = temp_contr[ij];
                ij++;
        }

        free(temp_contr);
        free(keyword);
        return contr;
#endif
    return NULL;
}

void Chkpt::wt_contr(double */*contr*/, const char */*key2*/)
{
#if 0
        int nprim;
        char *keyword;
        keyword = build_keyword("Contraction coefficients", key2);

        nprim = rd_nprim(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) contr,
                MAXANGMOM*nprim*sizeof(double));

        free(keyword);
#endif
}

extern "C" {
/*!
** chkpt_rd_contr(): Reads in the normalized contraction coefficients.
**
** Parameters: none
**
** Returns: double *contr Normalized contraction coefficients are
** returned as an array of doubles. In the checkpoint file they are
** stored as a matrix MAXANGMOM by the total number of primitives
** nprim, but each primitive Gaussian contributes to only one shell
** (and one basis function, of course), so most of these values are
** zero and not returned.
** \ingroup CHKPT
*/
        double *chkpt_rd_contr(void)
        {
                double *contr;
                contr = _default_chkpt_lib_->rd_contr();
                return contr;
        }

/*!
** chkpt_wt_contr(): Write out the normalized contraction coefficients.
**
** \param contr = The array of contraction coefficients.  The
**                ordering is that given in cints.
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_contr(double *contr, const char *key2)
        {
                _default_chkpt_lib_->wt_contr(contr, key2);
        }
}

