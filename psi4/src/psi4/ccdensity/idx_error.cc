/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {
namespace ccdensity {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
               int pq_sym, int rs_sym, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
  printer->Printf( "\n\tDPD Parameter Error in %s\n", message);
  printer->Printf(
          "\t-------------------------------------------------\n");
  printer->Printf(
          "\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
  printer->Printf("\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
          pq,rs,pq_sym,rs_sym);
  throw PsiException("ccdensity: error", __FILE__, __LINE__);
}


}} // namespace psi::ccdensity
