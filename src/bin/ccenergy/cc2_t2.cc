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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void dijabT2(void);
void cc2_faeT2(void);
void cc2_fmiT2(void);
void cc2_WmbijT2(void);
void cc2_WabeiT2(void);
void DT2(void);
void status(const char *s, std::string out);

void cc2_t2_build(void)
{

  DT2();

  if((params.ref == 0) || params.t2_coupled) { /** RHF or ROHF with coupled T2's **/ 
#ifdef TIME_CCENERGY
    timer_on("fT2", "outfile");
#endif
    cc2_faeT2(); cc2_fmiT2();
    if(params.print & 2) status("f -> T2", "outfile");
#ifdef TIME_CCENERGY
    timer_off("fT2", "outfile");
#endif
  }

#ifdef TIME_CCENERGY
  timer_on("WmbijT2", "outfile");
#endif
  cc2_WmbijT2();
  if(params.print & 2) status("Wmbij -> T2", "outfile");
#ifdef TIME_CCENERGY
  timer_off("WmbijT2", "outfile");
#endif

#ifdef TIME_CCENERGY
  timer_on("WabeiT2", "outfile");
#endif
  cc2_WabeiT2();
  if(params.print & 2) status("Wabei -> T2", "outfile");
#ifdef TIME_CCENERGY
  timer_off("WabeiT2", "outfile");
#endif

}
}} // namespace psi::ccenergy
