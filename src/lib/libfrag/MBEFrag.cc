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

#include "MBEFrag.h"
#include "FragOptions.h"
#include "Fragmenter.h"

namespace psi {
namespace LibFrag {

void MBEFrag::Copy(const MBEFrag& other) {
   this->Parents_=other.Parents_;
   this->MBEOrder_=other.MBEOrder_;
   this->Mult_=other.Mult_;
   this->Atoms_=other.Atoms_;
   this->Caps_=other.Caps_;
   this->Charges_=other.Charges_;
   this->Ghosts_=other.Ghosts_;
}

}
} //End namespaces

