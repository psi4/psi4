/*! \defgroup CHKPT libchkpt: The Checkpoint Interface library */

/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

/*!
** int *Chkpt::rd_am2canon_shell_order()
**
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** Returns: int *am2can_shell_order
** \ingroup CHKPT
*/
int *Chkpt::rd_am2canon_shell_order(const char *key2)
{
        int *am2can_sh_ord, nshell;
        char *keyword;
        keyword = build_keyword("AM -> canonical shell map", key2);

        nshell = rd_nshell(key2);
        am2can_sh_ord = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
                nshell*sizeof(int));

        free(keyword);
        return am2can_sh_ord;
}

/*!
** void Chkpt::wt_am2canon_shell_order()
**
** Writes out the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** \param am2can_shell_order = array to store the mapping array
**
** Returns: none
** \ingroup CHKPT
*/
void Chkpt::wt_am2canon_shell_order(int *am2can_sh_ord, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("AM -> canonical shell map", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
                nshell*sizeof(int));

  free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_am2canon_shell_order()
**
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** Returns: int *am2can_shell_order
** \ingroup CHKPT
*/
        int *chkpt_rd_am2canon_shell_order(void)
        {
                return _default_chkpt_lib_->rd_am2canon_shell_order();
        }

/*!
** void chkpt_wt_am2canon_shell_order()
**
** Writes out the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** \param am2can_shell_order = array to store the mapping array
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_am2canon_shell_order(int *am2can_sh_ord, const char *key2)
        {
                _default_chkpt_lib_->wt_am2canon_shell_order(am2can_sh_ord,
                  key2);
        }
}

