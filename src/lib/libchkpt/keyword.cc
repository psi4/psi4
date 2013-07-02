/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

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

