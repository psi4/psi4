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

// #include <libmoinfo/libmoinfo.h>
// 
// #include "blas.h"
// #include "index.h"
// #include "index_iterator.h"
// 
// namespace psi{ namespace psimrcc{
// 
// using namespace std;
// 
// int CCMatrixIterator::nirreps_=-1;
// 
// CCMatrixIterator::CCMatrixIterator(string left_str,string right_str,int irrep)
// {
//   ccindex_left  = blas->get_index(left_str);
//   ccindex_right = blas->get_index(right_str);
//   if(nirreps_<0) nirreps_ = moinfo->get_nirreps();
// /*  nelements_ = ccindex_->get_nelements();
//   element_irrep = ccindex_->get_element_irrep();
//   tuples     = ccindex_->get_tuples();
//   ind       = new short[nelements_];
//   ind_sym   = new int[nelements_];*/
//   reset();
// }
// 
// CCMatrixIterator::~CCMatrixIterator()
// {
//   delete[] ind;
//   delete[] ind_sym;
// }
// 
// void CCMatrixIterator::reset()
// {
//   sym = -1;
//   rel = -1;
//   abs = -1;
//   for(int n = 0; n < nelements_; ++n)
//     ind[n] = 0;
//   max_index_in_irrep_ = 0;
// }
// 
// int CCMatrixIterator::next_non_empty_irrep(int n)
// {
//   for(int h = n + 1; h < nirreps_; ++h)
//     if (ccindex_->get_pairpi(h) > 0) return h;
//   return -1;
// }
// 
// }}