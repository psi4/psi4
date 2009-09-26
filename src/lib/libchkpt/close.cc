/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

Chkpt::~Chkpt()
{
	psio->close(PSIF_CHKPT, 1);
	psio = NULL;
}

/*!
**  chkpt_close()  closes up the checkpoint file.
** 
**  Parameters: none, but chkpt_init must already have been called for 
**    this to work.  
**
**  Returns: none
**  \ingroup CHKPT
*/
extern "C" {
	int chkpt_close(void)
	{
		if (_default_chkpt_lib_) {
			delete _default_chkpt_lib_;
			_default_chkpt_lib_ = 0;
		}
		return 0;
	}
}
