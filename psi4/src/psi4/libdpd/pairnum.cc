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

#include "dpd.h"
using std::string;
using std::vector;
namespace psi {

int DPD::pairnum(string pair)
{
  vector<string> v = dpd_split(pair);

  int left, right;

  if(v.size() == 2) { // "pq"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[1] == moSpaces[i]) right = i;
    }
    if(left == right) return left*5; // unrestricted diagonal pair
    else if(left < right) return moSpaces.size()*5 + 2*(left*moSpaces.size()-left*(left+1)/2) + 2*(right-left-1);
    else if(left > right) return moSpaces.size()*5 + 2*(right*moSpaces.size()-right*(right+1)/2) + 2*(left-right-1) + 1;
  }
  else if(v.size() == 4) { // "p>q+" or "p>q-"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[2] == moSpaces[i]) right = i;
    }
    if(left != right) { throw; }
    if(v[3] == "+") return left*5 + 1;
    else if(v[3] == "-") return left*5 + 2;
  }
  else if(v.size() == 5) { // "p>=q+" or "p>=q-"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[3] == moSpaces[i]) right = i;
    }
    if(left != right) { throw; }
    if(v[4] == "+") return left*5 + 3;
    else if(v[4] == "-") return left*5 + 4;
  }
  return -1;
}

} // namespace psi
