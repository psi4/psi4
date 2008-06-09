/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

namespace psi {

char *Chkpt::build_keyword(char *key)
{
	char *keyword;
	int keylen;

	keylen = strlen(key) + strlen(chkpt_prefix) + 2;
	if(keylen > PSIO_KEYLEN) {
		printf("LIBCHKPT: requested key exceeds allowed LIBPSIO length: :%s:%s\n", 
			chkpt_prefix, key);
		exit(PSI_RETURN_FAILURE);
	}

	keyword = (char *) malloc((keylen+1)*sizeof(char));
	sprintf(keyword, ":%s:%s", chkpt_prefix, key);
	keyword[keylen] = '\0';

	return keyword;
}

	char *chkpt_build_keyword(char *key)
	{
		char *keyword;
		keyword = _default_chkpt_lib_->build_keyword(key);
		return keyword;
	}
}

