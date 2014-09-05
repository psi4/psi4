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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "psi4-dec.h"
#include "libparallel/ParallelPrinter.h"
namespace psi { namespace cclambda {

void status(const char *s, std::string out)
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
           boost::shared_ptr<OutFile>(new OutFile(out)));
  printer->Printf( "     %-15s...complete\n", s);

}

}} // namespace psi::cclambda
