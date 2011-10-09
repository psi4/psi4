/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

char *Chkpt::build_keyword(const char *key, const char *key2)
{
        char *keyword;
        int keylen;

        keylen = strlen(key) + strlen(key2) + strlen(chkpt_prefix) + 3;
        if(keylen > PSIO_KEYLEN) {
                printf("LIBCHKPT: requested key exceeds allowed LIBPSIO length: :%s:%s\n",
                        chkpt_prefix, key);
                exit(PSI_RETURN_FAILURE);
        }

        keyword = new char[keylen+1];
        if (key2[0] != '\0') {
          sprintf(keyword, ":%s:%s %s", chkpt_prefix, key, key2);
        }
        else {
          sprintf(keyword, ":%s:%s", chkpt_prefix, key);
        }
        keyword[keylen] = '\0';

        return keyword;
}

extern "C" {
        char *chkpt_build_keyword(const char *key)
        {
                char *keyword;
                keyword = _default_chkpt_lib_->build_keyword(key);
                return keyword;
        }
}

