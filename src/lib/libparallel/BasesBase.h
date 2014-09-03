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
#ifndef STUPIDWORKAROUND_H_
#define STUPIDWORKAROUND_H_

namespace psi{
///I needed to syphon off the calls that involve Worldcomm, hence this class
class BasesBase{
protected:
      ///Returns the integer of the Lucky MPI process
      int WhoIsSpecial()const;

      ///Returns true if this is the lucky MPI process that gets to read/write
      bool ImSpecial()const;
};

}


#endif /* STUPIDWORKAROUND_H_ */
