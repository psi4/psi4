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

const char *Chkpt::rd_label()
{
        const char *label;
        char *keyword;
        keyword = build_keyword("Label");

        label = (char *) malloc(80 * sizeof(char));

        psio->read_entry(PSIF_CHKPT, keyword, (char *) label, 80*sizeof(char));

        free(keyword);
        return label;
}

void Chkpt::wt_label(const char *label)
{
        char *keyword;
        keyword = build_keyword("Label");

        psio->write_entry(PSIF_CHKPT, keyword, (char*)label, 80*sizeof(char));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_label():  Reads the main chkpt label.
**
**   takes no arguments.
**
**   returns: pointer to the checkpoint label
** \ingroup CHKPT
*/
        const char *chkpt_rd_label(void)
        {
                const char *label;
                label = _default_chkpt_lib_->rd_label();
                return label;
        }

/*!
** chkpt_wt_label():  Writes the main chkpt label.
**
**  arguments:
**  \param label = The calculation label.
**
**   returns: none
** \ingroup CHKPT
*/

        void chkpt_wt_label(const char *label)
        {
                _default_chkpt_lib_->wt_label(label);
        }
}
