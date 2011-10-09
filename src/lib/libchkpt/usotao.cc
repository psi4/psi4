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

double **Chkpt::rd_usotao(const char *key2)
{
        double **usotao;
        int num_ao, num_so, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("SO->AO transmat", key2);

        num_ao = rd_nao(key2);
        num_so = rd_nso(key2);

        usotao = matrix<double>(num_so,num_ao);
        ptr = PSIO_ZERO;

        for(i=0;i<num_so;i++)
                psio->read(PSIF_CHKPT, keyword, (char *) usotao[i], (int) num_ao*sizeof(double), ptr, &ptr);

        free(keyword);
        return usotao;
}

void Chkpt::wt_usotao(double **usotao, const char *key2)
{
        int num_ao, num_so, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("SO->AO transmat", key2);

        num_ao = rd_nao(key2);
        num_so = rd_nso(key2);

        ptr = PSIO_ZERO;
        for(i=0;i<num_so;i++)
                psio->write(PSIF_CHKPT, keyword, (char *) usotao[i], (int) num_ao*sizeof(double), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_usotao(): Read in the SO to AO transformation matrix
**
** takes no arguments.
**
** returns: usotao = A num_so by num_ao matrix of doubles
**
** \ingroup CHKPT
*/
        double **chkpt_rd_usotao(void)
        {
                return _default_chkpt_lib_->rd_usotao();
        }

/*!
** chkpt_wt_usotao(): Writes out the SO to AO transformation matrix
**
** \param usotao = A num_so by num_ao matrix of doubles
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_usotao(double **usotao, const char *key2)
        {
                _default_chkpt_lib_->wt_usotao(usotao, key2);
        }
}
