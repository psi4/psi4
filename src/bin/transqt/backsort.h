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
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#ifndef _psi3_bin_transqt_backsort_h_
#define _psi3_bin_transqt_backsort_h_

namespace psi { namespace transqt {

void backsort_prep(int uhf);
void backsort(int first_tmp_file, double tolerance, int uhf);
void backsort_write(int i, int j, double **A, int kfirst, int klast,
		    int lfirst, int llast, int printflag, std::string OutFileRMR,
		       struct iwlbuf *twopdm_out, int uhf);

}} // end namespace psi::transqt
#endif // header guard
