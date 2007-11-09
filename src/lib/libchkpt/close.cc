/*!
  \file close.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
#include <libchkpt/chkpt.h>
}
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
**  arguments: none, but chkpt_init must already have been called for 
**    this to work.  
**
**  returns: zero.  Perhaps this, too, will change one day.
**  \ingroup (CHKPT)
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