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

char *Chkpt::rd_sym_label()
{
        char *sym_label;
        char *keyword;
        keyword = build_keyword("Symmetry label");

        sym_label = (char *) malloc(4*sizeof(char));

        psio->read_entry(PSIF_CHKPT, keyword, (char *) sym_label, 4*sizeof(char));

        sym_label[3] = '\0';

        free(keyword);
        return sym_label;
}

void Chkpt::wt_sym_label(char *sym_label)
{
        char *keyword;
        keyword = build_keyword("Symmetry label");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) sym_label, 4*sizeof(char));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_sym_label():  Reads in the symmetry label.
**
**   takes no arguments.
**
**   returns: symmetry = symmetry label.
**
** \ingroup CHKPT
*/
        char *chkpt_rd_sym_label(void)
        {
                return _default_chkpt_lib_->rd_sym_label();
        }

/*!
** chkpt_wt_sym_label():  Writes out the symmetry label.
**
** \param symmetry = symmetry label.
**
** returns none
**
** \ingroup CHKPT
*/
        void chkpt_wt_sym_label(char *sym_label)
        {
                _default_chkpt_lib_->wt_sym_label(sym_label);
        }
}

