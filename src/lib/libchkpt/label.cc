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
