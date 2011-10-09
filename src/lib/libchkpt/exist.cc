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

int Chkpt::exist(const char *keyword)
{
        int exists=0;

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL)
                exists=1;

        return exists;
}

extern "C" {
/*!
** chkpt_exist(): Checks to see if entry already exists in chkpt file. Note
** this function should be called only by functions in the chkpt library, as
** the calling function prepends the prefix.
**
**   takes no arguments.
**
**   returns: 1 if entry exists, 0 otherwise
**
** \ingroup CHKPT
*/
        int chkpt_exist(const char *keyword)
        {
                return(_default_chkpt_lib_->exist(keyword));
        }
}
