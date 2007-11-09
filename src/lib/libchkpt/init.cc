/*!
  \file init.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
#include <libchkpt/chkpt.h>
}
#include <libchkpt/chkpt.hpp>

using namespace psi;

/* Definition of global data */
Chkpt* psi::_default_chkpt_lib_ = 0;

//extern "C" {
/* first definition of chkpt_prefix */
//char chkpt_prefix[CHKPT_PREFIX_LEN];
//};

Chkpt::Chkpt(psi::PSIO *psioObject, int status) : psio(psioObject)
{
	char *prefix;
	psio_tocentry *this_entry;
	
	psio->open(PSIF_CHKPT, status);

	if(psio->tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
		prefix = rd_prefix();
		set_prefix(prefix);
		free(prefix);
	}
	else {
		set_prefix("");
		commit_prefix();  /* we assume that no default prefix existed in PSIF_CHKPT */
	}
}

void
Chkpt::rehash() {
  psio->rehash(PSIF_CHKPT);
}

extern "C" {
/*!
**  chkpt_init()  Initializes the checkpoint file for other chkpt_
**    functions to perform their duties.
**
**  arguments: 
**    int status: boolean indicating if the chkpt file should be
**                initialized (PSIO_OPEN_NEW) or the old chkpt 
**                file should be used (PSIO_OPEN_OLD).
**
**  returns: zero.  Perhaps this will change some day.
** \ingroup (CHKPT)
*/
	int chkpt_init(int status)
	{
		if (!_default_chkpt_lib_) {
			_default_chkpt_lib_ = new Chkpt(_default_psio_lib_, status);
			if (_default_chkpt_lib_ == 0) {
				fprintf(stderr, "LIBCHKPT::init() -- failed to allocate memory\n");
				exit(PSI_RETURN_FAILURE);
			}
		}
		return 0;
	}
}
