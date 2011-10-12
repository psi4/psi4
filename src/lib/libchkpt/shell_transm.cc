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

int **Chkpt::rd_shell_transm(const char *key2)
{
        int i, nshell, nirreps;
        int **shell_transm;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Shell transmat", key2);

        nshell = rd_nshell(key2);
        nirreps = rd_nirreps();

        shell_transm = matrix<int>(nshell,nirreps);
        ptr = PSIO_ZERO;
        for(i=0; i < nshell; i++)
                psio->read(PSIF_CHKPT, keyword, (char *) shell_transm[i],
                nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
        return shell_transm;
}

void Chkpt::wt_shell_transm(int **shell_transm, const char *key2)
{
        int i, nshell, nirreps;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Shell transmat", key2);

        nshell = rd_nshell(key2);
        nirreps = rd_nirreps();

        ptr = PSIO_ZERO;
        for(i=0; i < nshell; i++) {
                psio->write(PSIF_CHKPT, keyword, (char *) shell_transm[i],
                        nirreps*sizeof(int), ptr, &ptr);
        }

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_shell_transm():	Read in a matrix of nshell*nirreps integers
**			        that contains symmetry information.
**
**  takes no arguments.
**
**  returns:
**    shell_transm = matrix of nshell*nirrpes ints w/ symmetry info
**
** \ingroup CHKPT
*/
        int **chkpt_rd_shell_transm(void)
        {
                return _default_chkpt_lib_->rd_shell_transm();
        }

/*!
** chkpt_wt_shell_transm():	Write out a matrix of nshell*nirreps integers
**			        that contains symmetry information.
**
** \param shell_transm = matrix of nshell*nirreps ints w/ symmetry info
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_shell_transm(int **shell_transm, const char *key2)
        {
                _default_chkpt_lib_->wt_shell_transm(shell_transm, key2);
        }
}
